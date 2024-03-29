# Define LA (linear approximation-based) implementations of water models.

"""
    variable_flow(
        wm::AbstractLAModel;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true
    )

Creates flow-related variables related to linear approximation- (LA-) based
optimization models. First, creates flow variables for all node-connecting
components (e.g., pipes) in the network at subnetwork (or time) index `nw`,
e.g., `q_pipe[a]` for `a` in `pipe`. Then, creates continuous convex
combination variables used to construct necessary linear approximations, e.g.,
`lambda_pipe[a]` for `a` in `pipe`, where each is bounded between zero and one.
Finally, creates binary convex combination variables used for piecewise-linear
modeling of the approximating constraints, e.g., `x_pipe[a]` for `a` in `pipe`.
"""
function variable_flow(
    wm::AbstractLAModel;
    nw::Int = nw_id_default,
    bounded::Bool = true,
    report::Bool = true,
)
    for name in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
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


"""
    constraint_pipe_head_loss(
        wm::AbstractLAModel,
        n::Int,
        a::Int,
        node_fr::Int,
        node_to::Int,
        exponent::Float64,
        L::Float64,
        r::Float64,
        q_max_reverse::Float64,
        q_min_forward::Float64
    )

Adds constraints that model frictional head loss across a pipe via a piecewise-
linear approximation of the nonlinear, nonconvex head loss relationship. Here,
`wm` is the WaterModels object; `n` is the subnetwork (or time) index that is
considered; `a` is the index of the pipe; `node_fr` is the index of the tail
node of the pipe; `node_to` is the index of the head node of the pipe;
`exponent` is the exponent on flow in the head loss function (i.e., 1.852 for
Hazen-Williams head loss and 2.0 for Darcy-Weisbach head loss); `L` is the
length of the pipe; `r` is the resistance per unit length of the pipe;
`q_max_reverse` is the _maximum_ (negative) amount of flow when flow is
traveling in the negative direction (which corresponds to the _minimum_
magnitude of flow when traveling in the negative direction); and
`q_min_forward` is the _minimum_ (positive) amount of flow when flow is
traveling in the positive (forward) direction. Note that `q_max_reverse` and
`q_min_forward` are unused in this formulation since it is not direction-based.
"""
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


"""
    constraint_on_off_des_pipe_head_loss(
        wm::AbstractLAModel,
        n::Int,
        a::Int,
        node_fr::Int,
        node_to::Int,
        exponent::Float64,
        L::Float64,
        r::Float64,
        q_max_reverse::Float64,
        q_min_forward::Float64
    )

Add constraints that model frictional head loss across a design pipe via a
piecewise-linear approximation of the nonlinear, nonconvex head loss
relationship. Here, `wm` is the WaterModels object; `n` is the subnetwork (or
time) index that is considered; `a` is the index of the design pipe; `node_fr`
is the index of the tail node of the design pipe; `node_to` is the index of the
head node of the design pipe; `exponent` is the exponent on flow in the head
loss function (i.e., 1.852 for Hazen-Williams head loss and 2.0 for Darcy-
Weisbach head loss); `L` is the length of the design pipe; `r` is the
resistance per unit length of the design pipe; `q_max_reverse` is the _maximum_
(negative) amount of flow when flow is traveling in the negative direction
(which corresponds to the _minimum_ magnitude of flow when traveling in the
negative direction); and `q_min_forward` is the _minimum_ (positive) amount of
flow when flow is traveling in the positive (forward) direction. Note
direction-based flow limits are currently unused in these constraints.
"""
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


"""
    constraint_on_off_pump_head_gain(
        wm::AbstractLAModel,
        n::Int,
        a::Int,
        node_fr::Int,
        node_to::Int,
        q_min_forward::Float64
    )

Adds constraints that model the pump's head gain, if operating, as a piecewise-
linear approximation of a nonlinear function of flow rate. If the pump is
inactive, the head gain is restricted to a value of zero. Here, `wm` is the
WaterModels object, `n` is the subnetwork (or time) index that is considered,
`a` is the index of the pump, `node_fr` is the index of the tail node of the
pump, `node_to` is the index of the head node of the pump, and `q_min_forward`
is the _minimum_ (positive) amount of flow when the pump is active. Head gain
is assumed to be nonnegative and is directed from `node_fr` to `node_to`.
"""
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


"""
    constraint_on_off_pump_power(
        wm::AbstractLAModel,
        n::Int,
        a::Int,
        q_min_forward::Float64
    )

Adds constraints that model the pump's power consumption, if operating, as a
piecewise-linear approximation of a (potentially) nonlinear function of a more
accurate model. If the pump is inactive, the power is restricted to a value of
zero. Here, `wm` is the WaterModels object, `n` is the subnetwork (or time)
index that is considered, `a` is the index of the pump, and `q_min_forward` is
the _minimum_ (positive) amount of flow when the pump is active.
"""
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