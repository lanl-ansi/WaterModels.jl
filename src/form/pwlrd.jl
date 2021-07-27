# Define common PWLRD (piecewise linear relaxation- and direction-based) implementations of water
# distribution network constraints, which use directed flow variables.


########################################## VARIABLES ##########################################


"Create flow-related variables common to all directed flow models for node-connecting components."
function variable_flow(wm::AbstractPWLRDModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    for name in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
        # Create directed flow (`qp` and `qn`) variables for each component.
        _variable_component_flow(wm, name; nw = nw, bounded = bounded, report = report)

        # Create directed flow binary direction variables (`y`) for each component.
        _variable_component_direction(wm, name; nw = nw, report = report)
    end

    for name in ["des_pipe", "pipe"]
        # Create directed head difference (`dhp` and `dhn`) variables for each component.
        _variable_component_head_difference(wm, name; nw = nw, bounded = bounded, report = report)
    end
    
    # Create variables involved in convex combination constraints for pumps.
    pump_breakpoints = Dict{Int, Vector{Float64}}(a =>
        pump["flow_breakpoints"] for (a, pump) in ref(wm, nw, :pump))

    var(wm, nw)[:lambda_pump] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pump), k in 1:length(pump_breakpoints[a])],
        base_name = "$(nw)_lambda", lower_bound = 0.0, upper_bound = 1.0,
        start = comp_start_value(ref(wm, nw, :pump, a), "lambda_start", k))

    var(wm, nw)[:x_pump] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pump), k in 1:length(pump_breakpoints[a]) - 1],
        base_name = "$(nw)_x", binary = true,
        start = comp_start_value(ref(wm, nw, :pump, a), "x_start"))

    # Create variables involved in convex combination constraints for pipes.
    pipe_breakpoints_p = Dict{Int, Vector{Float64}}(a =>
        get_pipe_flow_breakpoints_positive(pipe)
        for (a, pipe) in ref(wm, nw, :pipe))

    var(wm, nw)[:lambda_p_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:length(pipe_breakpoints_p[a])],
        base_name = "$(nw)_lambda_p", lower_bound = 0.0, upper_bound = 1.0,
        start = comp_start_value(ref(wm, nw, :pipe, a), "lambda_p_start", k))

    var(wm, nw)[:x_p_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:length(pipe_breakpoints_p[a])-1],
        base_name = "$(nw)_x_p", binary = true,
        start = comp_start_value(ref(wm, nw, :pipe, a), "x_p_start"))

    pipe_breakpoints_n = Dict{Int, Vector{Float64}}(a =>
        sort(-get_pipe_flow_breakpoints_negative(pipe))
        for (a, pipe) in ref(wm, nw, :pipe))

    var(wm, nw)[:lambda_n_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:length(pipe_breakpoints_n[a])],
        base_name = "$(nw)_lambda_n", lower_bound = 0.0, upper_bound = 1.0,
        start = comp_start_value(ref(wm, nw, :pipe, a), "lambda_n_start", k))

    var(wm, nw)[:x_n_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:length(pipe_breakpoints_n[a])-1],
        base_name = "$(nw)_x_n", binary = true,
        start = comp_start_value(ref(wm, nw, :pipe, a), "x_n_start"))

    # Create variables involved in convex combination constraints for design pipes.
    des_pipe_breakpoints_p = Dict{Int, Vector{Float64}}(a =>
        get_pipe_flow_breakpoints_positive(des_pipe)
        for (a, des_pipe) in ref(wm, nw, :des_pipe))

    var(wm, nw)[:lambda_p_des_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :des_pipe), k in 1:length(des_pipe_breakpoints_p[a])],
        base_name = "$(nw)_lambda_p", lower_bound = 0.0, upper_bound = 1.0,
        start = comp_start_value(ref(wm, nw, :des_pipe, a), "lambda_p_start", k))

    var(wm, nw)[:x_p_des_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :des_pipe), k in 1:length(des_pipe_breakpoints_p[a])-1],
        base_name = "$(nw)_x_p", binary = true,
        start = comp_start_value(ref(wm, nw, :des_pipe, a), "x_p_start"))

    des_pipe_breakpoints_n = Dict{Int, Vector{Float64}}(a =>
        sort(-get_pipe_flow_breakpoints_negative(des_pipe))
        for (a, des_pipe) in ref(wm, nw, :des_pipe))

    var(wm, nw)[:lambda_n_des_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :des_pipe), k in 1:length(des_pipe_breakpoints_n[a])],
        base_name = "$(nw)_lambda_n", lower_bound = 0.0, upper_bound = 1.0,
        start = comp_start_value(ref(wm, nw, :des_pipe, a), "lambda_n_start", k))
    
    var(wm, nw)[:x_n_des_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :des_pipe), k in 1:length(des_pipe_breakpoints_n[a])-1],
        base_name = "$(nw)_x_n", binary = true,
        start = comp_start_value(ref(wm, nw, :des_pipe, a), "x_n_start"))
end


########################################## PIPES ##########################################


function constraint_pipe_head_loss(
    wm::AbstractPWLRDModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Get the variable for flow directionality.
    y = var(wm, n, :y_pipe, a)

    # Get variables for positive flow and head difference.
    qp, dhp = var(wm, n, :qp_pipe, a), var(wm, n, :dhp_pipe, a)

    # Get positively-directed convex combination variables.
    lambda_p, x_p = var(wm, n, :lambda_p_pipe), var(wm, n, :x_p_pipe)

    # Get the corresponding set of positive flow breakpoints.
    breakpoints_p = get_pipe_flow_breakpoints_positive(ref(wm, n, :pipe, a))
    bp_range, bp_range_m1 = 1:length(breakpoints_p), 1:length(breakpoints_p)-1

    # Add constraints for the positive flow piecewise approximation.
    c_1 = JuMP.@constraint(wm.model, sum(lambda_p[a, k] for k in bp_range) == y)
    qp_sum = sum(breakpoints_p[k] * lambda_p[a, k] for k in bp_range)
    c_2 = JuMP.@constraint(wm.model, qp_sum == qp)
    append!(con(wm, n, :pipe_head_loss, a), [c_1, c_2])

    if length(breakpoints_p) > 1
        # If there are multiple breakpoints, constrain the convex combination.
        c_3 = JuMP.@constraint(wm.model, sum(x_p[a, k] for k in bp_range_m1) == y)
        c_4 = JuMP.@constraint(wm.model, lambda_p[a, 1] <= x_p[a, 1])
        c_5 = JuMP.@constraint(wm.model, lambda_p[a, bp_range[end]] <= x_p[a, bp_range_m1[end]])
        append!(con(wm, n, :pipe_head_loss, a), [c_3, c_4, c_5])
    end

    # Add a constraint that upper-bounds the head loss variable.
    f_p = r .* breakpoints_p.^exponent
    f_p_ub_expr = sum(f_p[k] * lambda_p[a, k] for k in bp_range)
    c_6 = JuMP.@constraint(wm.model, dhp / L <= f_p_ub_expr)
    append!(con(wm, n, :pipe_head_loss, a), [c_6])

    # Loop over consequential points (i.e., those that have nonzero head loss).
    for flow_breakpoint in filter(x -> x > 0.0, breakpoints_p)
        # Add head loss outer (i.e., lower) approximations.
        lhs_p = _calc_head_loss_oa(qp, y, flow_breakpoint, exponent)
        c_7 = JuMP.@constraint(wm.model, r * lhs_p <= dhp / L)
        append!(con(wm, n, :pipe_head_loss, a), [c_7])
    end

    for k in 2:length(breakpoints_p)-1
        # Add the adjacency constraints for piecewise variables.
        c_8 = JuMP.@constraint(wm.model, lambda_p[a, k] <= x_p[a, k-1] + x_p[a, k])
        append!(con(wm, n, :pipe_head_loss, a), [c_8])
    end

    # Get variables for negative flow and head difference.
    qn, dhn = var(wm, n, :qn_pipe, a), var(wm, n, :dhn_pipe, a)

    # Get negatively-directed convex combination variable.
    lambda_n, x_n = var(wm, n, :lambda_n_pipe), var(wm, n, :x_n_pipe)

    # Get the corresponding set of negative flow breakpoints (negated).
    breakpoints_n = sort(-get_pipe_flow_breakpoints_negative(ref(wm, n, :pipe, a)))
    bn_range, bn_range_m1 = 1:length(breakpoints_n), 1:length(breakpoints_n)-1

    # Add constraints for the negative flow piecewise approximation.
    c_9 = JuMP.@constraint(wm.model, sum(lambda_n[a, k] for k in bn_range) == 1.0 - y)
    qn_sum = sum(breakpoints_n[k] * lambda_n[a, k] for k in bn_range)
    c_10 = JuMP.@constraint(wm.model, qn_sum == qn)
    append!(con(wm, n, :pipe_head_loss, a), [c_9, c_10])

    if length(breakpoints_n) > 1
        # If there are multiple breakpoints, constrain the convex combination.
        c_11 = JuMP.@constraint(wm.model, sum(x_n[a, k] for k in bn_range_m1) == 1.0 - y)
        c_12 = JuMP.@constraint(wm.model, lambda_n[a, 1] <= x_n[a, 1])
        c_13 = JuMP.@constraint(wm.model, lambda_n[a, bn_range[end]] <= x_n[a, bn_range_m1[end]])
        append!(con(wm, n, :pipe_head_loss, a), [c_11, c_12, c_13])
    end
    
    # Add a constraint that upper-bounds the head loss variable.
    f_n = r .* breakpoints_n.^exponent
    f_n_ub_expr = sum(f_n[k] * lambda_n[a, k] for k in bn_range)
    c_14 = JuMP.@constraint(wm.model, dhn / L <= f_n_ub_expr)
    append!(con(wm, n, :pipe_head_loss, a), [c_14])

    # Loop over consequential points (i.e., those that have nonzero head loss).
    for flow_breakpoint in filter(x -> x > 0.0, breakpoints_n)
        # Add head loss outer (i.e., lower) approximations.
        lhs_n = _calc_head_loss_oa(qn, 1.0 - y, flow_breakpoint, exponent)
        c_15 = JuMP.@constraint(wm.model, r * lhs_n <= dhn / L)
        append!(con(wm, n, :pipe_head_loss, a), [c_15])
    end

    for k in 2:length(breakpoints_n)-1
        # Add the adjacency constraints for piecewise variables.
        c_16 = JuMP.@constraint(wm.model, lambda_n[a, k] <= x_n[a, k-1] + x_n[a, k])
        append!(con(wm, n, :pipe_head_loss, a), [c_16])
    end
end


function constraint_on_off_des_pipe_head_loss(
    wm::AbstractPWLRDModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Gather design pipe flow direction and indicator variables.
    y, z = var(wm, n, :y_des_pipe, a), var(wm, n, :z_des_pipe, a)

    # Get variables for positive flow and head difference.
    qp, dhp = var(wm, n, :qp_des_pipe, a), var(wm, n, :dhp_des_pipe, a)

    # Gather positive convex combination variables.
    lambda_p, x_p = var(wm, n, :lambda_p_des_pipe), var(wm, n, :x_p_des_pipe)

    # Get the corresponding set of positive flow breakpoints.
    breakpoints_p = get_pipe_flow_breakpoints_positive(ref(wm, n, :des_pipe, a))
    bp_range, bp_range_m1 = 1:length(breakpoints_p), 1:length(breakpoints_p)-1

    # Add constraints for the positive flow piecewise approximation.
    c_1 = JuMP.@constraint(wm.model, sum(lambda_p[a, k] for k in bp_range) <= y)
    qp_sum = sum(breakpoints_p[k] * lambda_p[a, k] for k in bp_range)
    c_2 = JuMP.@constraint(wm.model, qp_sum == qp)
    append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_1, c_2])

    if length(breakpoints_p) > 1
        # If there are multiple breakpoints, constrain the convex combination.
        c_3 = JuMP.@constraint(wm.model, sum(x_p[a, k] for k in bp_range_m1) <= y)
        c_4 = JuMP.@constraint(wm.model, lambda_p[a, 1] <= x_p[a, 1])
        c_5 = JuMP.@constraint(wm.model, lambda_p[a, bp_range[end]] <= x_p[a, bp_range_m1[end]])
        append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_3, c_4, c_5])
    end

    # Add a constraint that upper-bounds the head loss variable.
    f_p = r .* breakpoints_p.^exponent
    f_p_ub_expr = sum(f_p[k] * lambda_p[a, k] for k in bp_range)
    c_6 = JuMP.@constraint(wm.model, dhp / L <= f_p_ub_expr)
    append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_6])

    # Loop over consequential points (i.e., those that have nonzero head loss).
    for flow_breakpoint in filter(x -> x > 0.0, breakpoints_p)
        # Add head loss outer (i.e., lower) approximations.
        lhs_p_1 = _calc_head_loss_oa(qp, y, flow_breakpoint, exponent)
        c_7 = JuMP.@constraint(wm.model, r * lhs_p_1 <= dhp / L)
        lhs_p_2 = _calc_head_loss_oa(qp, z, flow_breakpoint, exponent)
        c_8 = JuMP.@constraint(wm.model, r * lhs_p_2 <= dhp / L)
        append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_7, c_8])
    end

    for k in 2:length(breakpoints_p)-1
        # Add the adjacency constraints for piecewise variables.
        c_9 = JuMP.@constraint(wm.model, lambda_p[a, k] <= x_p[a, k-1] + x_p[a, k])
        append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_9])
    end

    # Get variables for negative flow and head difference.
    qn, dhn = var(wm, n, :qn_des_pipe, a), var(wm, n, :dhn_des_pipe, a)

    # Gather negative convex combination variables.
    lambda_n, x_n = var(wm, n, :lambda_n_des_pipe), var(wm, n, :x_n_des_pipe)

    # Get the corresponding set of negative flow breakpoints.
    breakpoints_n = sort(-get_pipe_flow_breakpoints_negative(ref(wm, n, :des_pipe, a)))
    bn_range, bn_range_m1 = 1:length(breakpoints_n), 1:length(breakpoints_n)-1

    # Add constraints for the negative flow piecewise approximation.
    c_10 = JuMP.@constraint(wm.model, sum(lambda_n[a, k] for k in bn_range) <= 1.0 - y)
    qn_sum = sum(breakpoints_n[k] * lambda_n[a, k] for k in bn_range)
    c_11 = JuMP.@constraint(wm.model, qn_sum == qn)
    append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_10, c_11])

    if length(breakpoints_n) > 1
        # If there are multiple breakpoints, constrain the convex combination.
        c_12 = JuMP.@constraint(wm.model, sum(x_n[a, k] for k in bn_range_m1) <= 1.0 - y)
        c_13 = JuMP.@constraint(wm.model, lambda_n[a, 1] <= x_n[a, 1])
        c_14 = JuMP.@constraint(wm.model, lambda_n[a, bn_range[end]] <= x_n[a, bn_range_m1[end]])
        append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_12, c_13, c_14])
    end

    # Add a constraint that upper-bounds the head loss variable.
    f_n = r .* breakpoints_n.^exponent
    f_n_ub_expr = sum(f_n[k] * lambda_n[a, k] for k in bn_range)
    c_15 = JuMP.@constraint(wm.model, dhn / L <= f_n_ub_expr)
    append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_15])

    # Loop over consequential points (i.e., those that have nonzero head loss).
    for flow_breakpoint in filter(x -> x > 0.0, breakpoints_n)
        # Add head loss outer (i.e., lower) approximations.
        lhs_n_1 = _calc_head_loss_oa(qn, 1.0 - y, flow_breakpoint, exponent)
        c_16 = JuMP.@constraint(wm.model, r * lhs_n_1 <= dhn / L)
        lhs_n_2 = _calc_head_loss_oa(qn, z, flow_breakpoint, exponent)
        c_17 = JuMP.@constraint(wm.model, r * lhs_n_2 <= dhn / L)
        append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_16, c_17])
    end

    for k in 2:length(breakpoints_n)-1
        # Add the adjacency constraints for piecewise variables.
        c_18 = JuMP.@constraint(wm.model, lambda_n[a, k] <= x_n[a, k-1] + x_n[a, k])
        append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_18])
    end

    # Constrain the negative and positive piecewise variables per status.
    x_p_sum = length(bp_range_m1) > 0 ? sum(x_p[a, k] for k in bp_range_m1) : 0.0
    x_n_sum = length(bn_range_m1) > 0 ? sum(x_n[a, k] for k in bn_range_m1) : 0.0
    c_19 = JuMP.@constraint(wm.model, x_p_sum + x_n_sum == z)

    lambda_p_sum = sum(lambda_p[a, k] for k in bp_range)
    lambda_n_sum = sum(lambda_n[a, k] for k in bn_range) 
    c_20 = JuMP.@constraint(wm.model, lambda_p_sum + lambda_n_sum == z)
    append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_19, c_20])
end


########################################## PUMPS ##########################################


"Add constraints associated with modeling a pump's head gain."
function constraint_on_off_pump_head_gain(wm::AbstractPWLRDModel, n::Int, a::Int, node_fr::Int, node_to::Int, q_min_forward::Float64)
    # Gather pump flow, head gain, and status variables.
    qp, g, z = var(wm, n, :qp_pump, a), var(wm, n, :g_pump, a), var(wm, n, :z_pump, a)

    # Gather convex combination variables.
    lambda, x = var(wm, n, :lambda_pump), var(wm, n, :x_pump)

    # Get the corresponding set of flow breakpoints.
    @assert haskey(ref(wm, n, :pump, a), "flow_breakpoints")
    breakpoints = ref(wm, n, :pump, a, "flow_breakpoints")
    bp_range, bp_range_m1 = 1:length(breakpoints), 1:length(breakpoints)-1

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, k] for k in bp_range) == z)
    c_2 = JuMP.@constraint(wm.model, sum(x[a, k] for k in bp_range_m1) == z)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, bp_range[end]] <= x[a, bp_range_m1[end]])

    # Add a constraint for the flow piecewise approximation.
    qp_lhs = sum(breakpoints[k] * lambda[a, k] for k in bp_range)
    c_5 = JuMP.@constraint(wm.model, qp_lhs == qp)

    # Get pump head curve function and its derivative.
    head_curve_function = ref(wm, n, :pump, a, "head_curve_function")
    head_curve_derivative = ref(wm, n, :pump, a, "head_curve_derivative")

    # Add a constraint that lower-bounds the head gain variable.
    f_all = head_curve_function.(breakpoints)
    gain_lb_expr = sum(f_all[k] .* lambda[a, k] for k in bp_range)
    c_6 = JuMP.@constraint(wm.model, gain_lb_expr <= g)

    # Append the constraint array.
    append!(con(wm, n, :on_off_pump_head_gain, a), [c_1, c_2, c_3, c_4, c_5, c_6])

    for flow_breakpoint in filter(q_val -> q_val > 0.0, breakpoints)
        # Add head gain outer (i.e., upper) approximations.
        f = head_curve_function(flow_breakpoint)
        df = head_curve_derivative(flow_breakpoint)
        f_rhs = f * z + df * (qp - flow_breakpoint * z)
        c = JuMP.@constraint(wm.model, g <= f_rhs)
        append!(con(wm, n, :on_off_pump_head_gain, a), [c])
    end

    for k in 2:length(breakpoints)-1
        # Add the adjacency constraints for piecewise variables.
        adjacency = x[a, k-1] + x[a, k]
        c = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :on_off_pump_head_gain, a), [c])
    end
end


function constraint_on_off_pump_power(wm::AbstractPWLRDModel, n::Int, a::Int, q_min_forward::Float64)
    # Gather pump power and status variables.
    P, lambda = var(wm, n, :P_pump, a), var(wm, n, :lambda_pump)

    # Get the corresponding set of flow breakpoints.
    @assert haskey(ref(wm, n, :pump, a), "flow_breakpoints")
    breakpoints = ref(wm, n, :pump, a, "flow_breakpoints")
    bp_range = 1:length(breakpoints)
    
    # Add a constraint that lower-bounds the power variable.
    f_ua = _calc_pump_power_ua(wm, n, a, breakpoints)
    power_lb_expr = sum(f_ua[k] * lambda[a, k] for k in bp_range)
    c_1 = JuMP.@constraint(wm.model, power_lb_expr <= P)

    # Add a constraint that upper-bounds the power variable.
    f_oa = _calc_pump_power_oa(wm, n, a, breakpoints)
    power_ub_expr = sum(f_oa[k] * lambda[a, k] for k in bp_range)
    c_2 = JuMP.@constraint(wm.model, P <= power_ub_expr)

    # Append the :on_off_pump_power constraint array.
    append!(con(wm, n, :on_off_pump_power)[a], [c_1, c_2])
end