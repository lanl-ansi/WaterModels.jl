# Define LA (linear approximation-based) implementations of water models.

"Creates flow variables for `LA` formulations (`q`, `lambda`, `x_pw`)."
function variable_flow(
    wm::AbstractLAModel;
    nw::Int = nw_id_default,
    bounded::Bool = true,
    report::Bool = true,
)
    for name in _LINK_COMPONENTS
        # Create flow variables for each node-connecting component.
        _variable_component_flow(wm, name; nw = nw, bounded = bounded, report = report)
    end

    # Create weights involved in convex combination constraints for pipes.
    var(wm, nw)[:lambda_pipe] = JuMP.@variable(
        wm.model,
        [
            a in ids(wm, nw, :pipe),
            k in 1:length(ref(wm, nw, :pipe, a, "flow_partition"))
        ],
        base_name = "$(nw)_lambda",
        lower_bound = 0.0,
        upper_bound = 1.0,
        start = comp_start_value(ref(wm, nw, :pipe, a), "lambda_start", k)
    )

    # Create weights involved in convex combination constraints for design pipes.
    var(wm, nw)[:lambda_des_pipe] = JuMP.@variable(
        wm.model,
        [
            a in ids(wm, nw, :des_pipe),
            k in 1:length(ref(wm, nw, :des_pipe, a, "flow_partition")),
        ],
        base_name = "$(nw)_lambda",
        lower_bound = 0.0,
        upper_bound = 1.0,
        start = comp_start_value(ref(wm, nw, :des_pipe, a), "lambda_start", k)
    )

    # Create weights involved in convex combination constraints for pumps.
    var(wm, nw)[:lambda_pump] = JuMP.@variable(
        wm.model,
        [
            a in ids(wm, nw, :pump),
            k in 1:length(ref(wm, nw, :pump, a, "flow_partition"))
        ],
        base_name = "$(nw)_lambda",
        lower_bound = 0.0,
        upper_bound = 1.0,
        start = comp_start_value(ref(wm, nw, :pump, a), "lambda_start", k)
    )

    # Create binary variables for pipe convex combination constraints.
    var(wm, nw)[:x_pipe] = JuMP.@variable(
        wm.model,
        [
            a in ids(wm, nw, :pipe),
            k in 1:length(ref(wm, nw, :pipe, a, "flow_partition"))-1,
        ],
        base_name = "$(nw)_x_pipe",
        binary = true,
        start = comp_start_value(ref(wm, nw, :pipe, a), "x_start")
    )

    # Create binary variables for design pipe convex combination constraints.
    var(wm, nw)[:x_des_pipe] = JuMP.@variable(
        wm.model,
        [
            a in ids(wm, nw, :des_pipe),
            k in 1:length(ref(wm, nw, :des_pipe, a, "flow_partition"))-1,
        ],
        base_name = "$(nw)_x_des_pipe",
        binary = true,
        start = comp_start_value(ref(wm, nw, :des_pipe, a), "x_start")
    )

    # Create binary variables for pump convex combination constraints.
    var(wm, nw)[:x_pump] = JuMP.@variable(
        wm.model,
        [
            a in ids(wm, nw, :pump),
            k in 1:length(ref(wm, nw, :pump, a, "flow_partition"))-1,
        ],
        base_name = "$(nw)_x_pump",
        binary = true,
        start = comp_start_value(ref(wm, nw, :pump, a), "x_start")
    )
end


function constraint_pipe_head_loss(
    wm::AbstractLAModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    exponent::Float64,
    L::Float64,
    r::Float64,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # Get required variables.
    q, h_i, h_j = var(wm, n, :q_pipe, a), var(wm, n, :h, node_fr), var(wm, n, :h, node_to)
    lambda, x = var(wm, n, :lambda_pipe), var(wm, n, :x_pipe)

    # Get the corresponding flow partitioning.
    @assert haskey(ref(wm, n, :pipe, a), "flow_partition")
    partition = ref(wm, n, :pipe, a, "flow_partition")
    bp_range, bp_range_m1 = 1:length(partition), 1:length(partition)-1

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, k] for k in bp_range) == 1.0)
    c_2 = JuMP.@constraint(wm.model, sum(x[a, k] for k in bp_range_m1) == 1.0)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, bp_range[end]] <= x[a, bp_range_m1[end]])

    # Add a constraint for the flow piecewise approximation.
    q_lhs = sum(partition[k] * lambda[a, k] for k in bp_range)
    c_5 = JuMP.@constraint(wm.model, q_lhs == q)

    # Add a constraint for the head loss piecewise approximation.
    f = _calc_head_loss_values(partition, exponent)
    lhs = r * sum(f[k] * lambda[a, k] for k in bp_range)
    c_6 = JuMP.@constraint(wm.model, lhs == (h_i - h_j) / L)

    # Append the constraint array with the above-generated constraints.
    append!(con(wm, n, :pipe_head_loss, a), [c_1, c_2, c_3, c_4, c_5, c_6])

    # Add the adjacency constraints.
    for k = 2:length(partition)-1
        adjacency = x[a, k-1] + x[a, k]
        c_7_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :pipe_head_loss, a), [c_7_k])
    end
end


function constraint_on_off_des_pipe_head_loss(
    wm::AbstractLAModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    exponent::Float64,
    L::Float64,
    r::Float64,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # Get required variables.
    q, z = var(wm, n, :q_des_pipe, a), var(wm, n, :z_des_pipe, a)
    dh = var(wm, n, :dh_des_pipe, a) # Zero when design pipe is not selected.
    lambda, x = var(wm, n, :lambda_des_pipe), var(wm, n, :x_des_pipe)

    # Get the corresponding flow partitioning.
    @assert haskey(ref(wm, n, :des_pipe, a), "flow_partition")
    partition = ref(wm, n, :des_pipe, a, "flow_partition")
    bp_range, bp_range_m1 = 1:length(partition), 1:length(partition)-1

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, k] for k in bp_range) == z)
    c_2 = JuMP.@constraint(wm.model, sum(x[a, k] for k in bp_range_m1) == z)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, bp_range[end]] <= x[a, bp_range_m1[end]])

    # Add a constraint for the flow piecewise approximation.
    q_lhs = sum(partition[k] * lambda[a, k] for k in bp_range)
    c_5 = JuMP.@constraint(wm.model, q_lhs == q)

    # Add a constraint for the head loss piecewise approximation.
    f = _calc_head_loss_values(partition, exponent)
    lhs = r * sum(f[k] * lambda[a, k] for k in bp_range)

    # TODO: Use a McCormick expansion of the below multiplication with z.
    c_6 = JuMP.@constraint(wm.model, lhs == dh / L)

    # Append the constraint array with the above-generated constraints.
    append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_1, c_2, c_3, c_4, c_5, c_6])

    # Add the adjacency constraints.
    for k = 2:length(partition)-1
        adjacency = x[a, k-1] + x[a, k]
        c_7_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_7_k])
    end
end


function constraint_on_off_pump_head_gain(
    wm::AbstractLAModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    q_min_forward::Float64,
)
    # Gather common variables.
    q, g, z = var(wm, n, :q_pump, a), var(wm, n, :g_pump, a), var(wm, n, :z_pump, a)
    lambda, x = var(wm, n, :lambda_pump), var(wm, n, :x_pump)

    # Get the corresponding flow partitioning.
    @assert haskey(ref(wm, n, :pump, a), "flow_partition")
    partition = ref(wm, n, :pump, a, "flow_partition")
    bp_range, bp_range_m1 = 1:length(partition), 1:length(partition)-1

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, k] for k in bp_range) == z)
    c_2 = JuMP.@constraint(wm.model, sum(x[a, k] for k in bp_range_m1) == z)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, bp_range[end]] <= x[a, bp_range_m1[end]])

    # Add a constraint for the flow piecewise approximation.
    q_lhs = sum(partition[k] * lambda[a, k] for k in bp_range)
    c_5 = JuMP.@constraint(wm.model, q_lhs == q)

    # Add a constraint for the head gain piecewise approximation.
    head_curve_function = _calc_head_curve_function(ref(wm, n, :pump, a))

    # Add a constraint that linearly approximates the head gain variable.
    f = head_curve_function.(partition)
    g_lhs = sum(f[k] * lambda[a, k] for k in bp_range)
    c_6 = JuMP.@constraint(wm.model, g_lhs == g)

    # Append the constraint array with the above-generated constraints.
    append!(con(wm, n, :on_off_pump_head_gain, a), [c_1, c_2, c_3, c_4, c_5, c_6])

    for k = 2:length(partition)-1
        # Add adjacency constraints for each interval.
        adjacency = x[a, k-1] + x[a, k]
        c_7_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :on_off_pump_head_gain, a), [c_7_k])
    end
end


function constraint_on_off_pump_power(
    wm::AbstractLAModel,
    n::Int,
    a::Int,
    q_min_forward::Float64,
)
    # Gather pump power and convex combination variables.
    P, lambda = var(wm, n, :P_pump, a), var(wm, n, :lambda_pump)

    # Get the corresponding flow partition.
    @assert haskey(ref(wm, n, :pump, a), "flow_partition")
    partition = ref(wm, n, :pump, a, "flow_partition")

    # Generate the set of power approximation points.
    f_all = _calc_pump_power(wm, n, a, partition)

    # Add a constraint that lower-bounds the power variable.
    power_expr = sum(f_all[k] * lambda[a, k] for k = 1:length(partition))
    c = JuMP.@constraint(wm.model, power_expr == P)

    # Append the :on_off_pump_power constraint array.
    append!(con(wm, n, :on_off_pump_power)[a], [c])
end