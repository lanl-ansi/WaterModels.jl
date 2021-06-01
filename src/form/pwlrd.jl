# Define common PWLRD (piecewise linear relaxation- and direction-based) implementations of water
# distribution network constraints, which use directed flow variables.


########################################## VARIABLES ##########################################


"Create flow-related variables common to all directed flow models for node-connecting components."
function variable_flow(wm::AbstractPWLRDModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    for name in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
        # Create directed flow (`qp` and `qn`) variables for each component.
        _variable_component_flow(wm, name; nw=nw, bounded=bounded, report=report)

        # Create directed flow binary direction variables (`y`) for each component.
        _variable_component_direction(wm, name; nw=nw, report=report)
    end

    for name in ["des_pipe", "pipe"]
        # Create directed head difference (`dhp` and `dhn`) variables for each component.
        _variable_component_head_difference(wm, name; nw=nw, bounded=bounded, report=report)
    end

    # Get the number of breakpoints for pumps.
    pump_breakpoints = Dict{Int, Vector{Float64}}(a =>
        get_pump_flow_lower_breakpoints_positive(pump)
        for (a, pump) in ref(wm, nw, :pump))
    
    # Create weights involved in convex combination constraints for pumps.
    var(wm, nw)[:lambda_pump] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pump), k in 1:length(pump_breakpoints[a])],
        base_name="$(nw)_lambda", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :pump, a), "lambda_start", k))

    # Create binary variables involved in convex combination constraints for pumps.
    var(wm, nw)[:x_pw_pump] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pump), k in 1:length(pump_breakpoints[a])-1], base_name="$(nw)_x_pw",
        binary=true, start=comp_start_value(ref(wm, nw, :pump, a), "x_pw_start"))

    # Create weights involved in convex combination constraints for pipes.
    pipe_breakpoints_p = Dict{Int, Vector{Float64}}(a =>
        get_pipe_flow_upper_breakpoints_positive(pipe)
        for (a, pipe) in ref(wm, nw, :pipe))

    var(wm, nw)[:lambda_p_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:length(pipe_breakpoints_p[a])],
        base_name="$(nw)_lambda_p", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :pipe, a), "lambda_p_start", k))

    # Create weights involved in convex combination constraints for pipes.
    pipe_breakpoints_n = Dict{Int, Vector{Float64}}(a =>
        get_pipe_flow_upper_breakpoints_negative(pipe)
        for (a, pipe) in ref(wm, nw, :pipe))

    var(wm, nw)[:lambda_n_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:length(pipe_breakpoints_n[a])],
        base_name="$(nw)_lambda_n", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :pipe, a), "lambda_n_start", k))

    # Create weights involved in convex combination constraints for design pipes.
    des_pipe_breakpoints_p = Dict{Int, Vector{Float64}}(a =>
        get_pipe_flow_upper_breakpoints_positive(des_pipe)
        for (a, des_pipe) in ref(wm, nw, :des_pipe))

    var(wm, nw)[:lambda_p_des_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :des_pipe), k in 1:length(des_pipe_breakpoints_p[a])],
        base_name="$(nw)_lambda_p", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :des_pipe, a), "lambda_p_start", k))

    # Create weights involved in convex combination constraints for design pipes.
    des_pipe_breakpoints_n = Dict{Int, Vector{Float64}}(a =>
        get_pipe_flow_upper_breakpoints_negative(des_pipe)
        for (a, des_pipe) in ref(wm, nw, :des_pipe))

    var(wm, nw)[:lambda_n_des_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :des_pipe), k in 1:length(des_pipe_breakpoints_n[a])],
        base_name="$(nw)_lambda_n", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :des_pipe, a), "lambda_n_start", k))
    
    # Create binary variables involved in convex combination constraints for pipes.
    var(wm, nw)[:x_p_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:length(pipe_breakpoints_p[a])-1], base_name="$(nw)_x_p",
        binary=true, start=comp_start_value(ref(wm, nw, :pipe, a), "x_p_start"))

    # Create binary variables involved in convex combination constraints for pipes.
    var(wm, nw)[:x_n_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:length(pipe_breakpoints_n[a])-1], base_name="$(nw)_x_n",
        binary=true, start=comp_start_value(ref(wm, nw, :pipe, a), "x_n_start"))

    # Create binary variables involved in convex combination constraints for design pipes.
    var(wm, nw)[:x_p_des_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :des_pipe), k in 1:length(des_pipe_breakpoints_p[a])-1], base_name="$(nw)_x_p",
        binary=true, start=comp_start_value(ref(wm, nw, :des_pipe, a), "x_p_start"))

    # Create binary variables involved in convex combination constraints for design pipes.
    var(wm, nw)[:x_n_des_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :des_pipe), k in 1:length(des_pipe_breakpoints_n[a])-1], base_name="$(nw)_x_n",
        binary=true, start=comp_start_value(ref(wm, nw, :des_pipe, a), "x_n_start"))
end


########################################## PIPES ##########################################


function constraint_pipe_head_loss(
    wm::AbstractPWLRDModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Get object from the WaterModels reference dictionary.
    pipe = ref(wm, n, :pipe, a)

    # Gather directed pipe flow and direction variables.
    qp, qn, y = var(wm, n, :qp_pipe, a), var(wm, n, :qn_pipe, a), var(wm, n, :y_pipe, a)

    # Gather directed head loss variables.
    dhp, dhn = var(wm, n, :dhp_pipe, a), var(wm, n, :dhn_pipe, a)

    # Add a constraint for the flow piecewise approximation.
    lambda_p, x_p = var(wm, n, :lambda_p_pipe), var(wm, n, :x_p_pipe)
    breakpoints_p = get_pipe_flow_upper_breakpoints_positive(pipe)

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda_p[a, k] for k in 1:length(breakpoints_p)) == y)
    qp_lhs = sum(qp_hat * lambda_p[a, k] for (k, qp_hat) in enumerate(breakpoints_p))
    c_2 = JuMP.@constraint(wm.model, qp_lhs == qp)
    append!(con(wm, n, :pipe_head_loss, a), [c_1, c_2])

    if length(breakpoints_p) > 1
        c_3 = JuMP.@constraint(wm.model, sum(x_p[a, k] for k in 1:length(breakpoints_p)-1) == y)
        c_4 = JuMP.@constraint(wm.model, lambda_p[a, 1] <= x_p[a, 1])
        c_5 = JuMP.@constraint(wm.model, lambda_p[a, length(breakpoints_p)] <= x_p[a, length(breakpoints_p)-1])
        append!(con(wm, n, :pipe_head_loss, a), [c_3, c_4, c_5])
    end

    # Add a constraint that upper-bounds the head loss variable.
    f_p = r .* breakpoints_p.^exponent
    loss_p_ub_expr = sum(f_p[k] * lambda_p[a, k] for k in 1:length(breakpoints_p))
    c_6 = JuMP.@constraint(wm.model, dhp / L <= loss_p_ub_expr)
    append!(con(wm, n, :pipe_head_loss, a), [c_6])

    for qp_hat in get_pipe_flow_lower_breakpoints_positive(pipe)
        # Add head loss outer (i.e., lower) approximations.
        lhs_p = _calc_head_loss_oa(qp, y, qp_hat, exponent)
        c = JuMP.@constraint(wm.model, r * lhs_p <= dhp / L)
        append!(con(wm, n, :pipe_head_loss, a), [c])
    end

    for k in 2:length(breakpoints_p)-1
        # Add the adjacency constraints for piecewise variables.
        c = JuMP.@constraint(wm.model, lambda_p[a, k] <= x_p[a, k-1] + x_p[a, k])
        append!(con(wm, n, :pipe_head_loss, a), [c])
    end

    # Add a constraint for the flow piecewise approximation.
    lambda_n, x_n = var(wm, n, :lambda_n_pipe), var(wm, n, :x_n_pipe)
    breakpoints_n = -get_pipe_flow_upper_breakpoints_negative(pipe)
    qn_lhs = sum(qn_hat * lambda_n[a, k] for (k, qn_hat) in enumerate(breakpoints_n))

    c_7 = JuMP.@constraint(wm.model, sum(lambda_n[a, k] for k in 1:length(breakpoints_n)) == 1.0 - y)
    c_8 = JuMP.@constraint(wm.model, qn_lhs == qn)
    append!(con(wm, n, :pipe_head_loss, a), [c_7, c_8])

    if length(breakpoints_n) > 1
        c_9 = JuMP.@constraint(wm.model, lambda_n[a, 1] <= x_n[a, 1])
        c_10 = JuMP.@constraint(wm.model, sum(x_n[a, k] for k in 1:length(breakpoints_n)-1) == 1.0 - y)
        c_11 = JuMP.@constraint(wm.model, lambda_n[a, length(breakpoints_n)] <= x_n[a, length(breakpoints_n)-1])
        append!(con(wm, n, :pipe_head_loss, a), [c_9, c_10, c_11])
    end
    
    # Add a constraint that upper-bounds the head loss variable.
    f_n =  r .* breakpoints_n.^exponent
    loss_n_ub_expr = sum(f_n[k] .* lambda_n[a, k] for (k, qn_hat) in enumerate(breakpoints_n))
    c_12 = JuMP.@constraint(wm.model, dhn / L <= loss_n_ub_expr)
    append!(con(wm, n, :pipe_head_loss, a), [c_12])

    for qn_hat in -get_pipe_flow_lower_breakpoints_negative(pipe)
        # Add head loss outer (i.e., lower) approximations.
        lhs_n = _calc_head_loss_oa(qn, 1.0 - y, qn_hat, exponent)
        c = JuMP.@constraint(wm.model, r * lhs_n <= dhn / L)
        append!(con(wm, n, :pipe_head_loss, a), [c])
    end

    for k in 2:length(breakpoints_n)-1
        # Add the adjacency constraints for piecewise variables.
        c = JuMP.@constraint(wm.model, lambda_n[a, k] <= x_n[a, k-1] + x_n[a, k])
        append!(con(wm, n, :pipe_head_loss, a), [c])
    end
end


function constraint_on_off_des_pipe_head_loss(
    wm::AbstractPWLRDModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Get the number of breakpoints for the pipe.
    pipe_breakpoints = get(wm.ext, :pipe_breakpoints, 2)

    # Gather directed pipe flow and direction variables.
    qp, qn = var(wm, n, :qp_des_pipe, a), var(wm, n, :qn_des_pipe, a)
    y, z = var(wm, n, :y_des_pipe, a), var(wm, n, :z_des_pipe, a)

    # Gather directed head loss variables.
    dhp, dhn = var(wm, n, :dhp_des_pipe, a), var(wm, n, :dhn_des_pipe, a)

    # Gather convex combination variables.
    lambda_p, x_p = var(wm, n, :lambda_p_des_pipe), var(wm, n, :x_p_des_pipe)
    lambda_n, x_n = var(wm, n, :lambda_n_des_pipe), var(wm, n, :x_n_des_pipe)

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda_p[a, :]) <= y)
    c_2 = JuMP.@constraint(wm.model, sum(x_p[a, :]) <= y)
    c_3 = JuMP.@constraint(wm.model, lambda_p[a, 1] <= x_p[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda_p[a, end] <= x_p[a, end])

    c_5 = JuMP.@constraint(wm.model, sum(lambda_n[a, :]) <= 1.0 - y)
    c_6 = JuMP.@constraint(wm.model, sum(x_n[a, :]) <= 1.0 - y)
    c_7 = JuMP.@constraint(wm.model, lambda_n[a, 1] <= x_n[a, 1])
    c_8 = JuMP.@constraint(wm.model, lambda_n[a, end] <= x_n[a, end])

    # Constrain the negative and positive piecewise variables per status.
    c_9 = JuMP.@constraint(wm.model, sum(x_p[a, :]) + sum(x_n[a, :]) == z)
    c_10 = JuMP.@constraint(wm.model, sum(lambda_p[a, :]) + sum(lambda_n[a, :]) == z)

    # Add a constraint for the flow piecewise approximation.
    breakpoints_p = range(max(0.0, q_min_forward), stop=JuMP.upper_bound(qp), length=pipe_breakpoints)
    qp_lhs = sum(breakpoints_p[k] * lambda_p[a, k] for k in 1:pipe_breakpoints)
    c_5 = JuMP.@constraint(wm.model, qp_lhs == qp)

    # Add a constraint that upper-bounds the head loss variable.
    f_p = r .* breakpoints_p.^exponent
    loss_p_ub_expr = sum(f_p[k] * lambda_p[a, k] for k in 1:pipe_breakpoints)
    c_6 = JuMP.@constraint(wm.model, inv(L) * dhp <= loss_p_ub_expr)

    # Add a constraint for the flow piecewise approximation.
    breakpoints_n = range(max(0.0, -q_max_reverse), stop=JuMP.upper_bound(qn), length=pipe_breakpoints)
    qn_lhs = sum(breakpoints_n[k] * lambda_n[a, k] for k in 1:pipe_breakpoints)
    c_5 = JuMP.@constraint(wm.model, qn_lhs == qn)

    # Add a constraint that upper-bounds the head loss variable.
    f_n =  r .* breakpoints_n.^exponent
    loss_n_ub_expr = sum(f_n[k] .* lambda_n[a, k] for k in 1:pipe_breakpoints)
    c_6 = JuMP.@constraint(wm.model, inv(L) * dhn <= loss_n_ub_expr)

    # Append the constraint array.
    append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_1, c_2, c_3, c_4, c_5, c_6])

    for qp_hat in breakpoints_p
        # Add head loss outer (i.e., lower) approximations.
        lhs_p_y = _calc_head_loss_oa(qp, y, qp_hat, exponent)
        c_7_k_y = JuMP.@constraint(wm.model, r * lhs_p_y <= dhp / L)

        lhs_p_z = _calc_head_loss_oa(qp, z, qp_hat, exponent)
        c_7_k_z = JuMP.@constraint(wm.model, r * lhs_p_z <= dhp / L)

        append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_7_k_y, c_7_k_z])
    end

    for qn_hat in breakpoints_n
        # Add head loss outer (i.e., lower) approximations.
        lhs_n_y = _calc_head_loss_oa(qn, 1.0 - y, qn_hat, exponent)
        c_7_k_y = JuMP.@constraint(wm.model, r * lhs_n_y <= dhn / L)

        lhs_n_z = _calc_head_loss_oa(qn, z, qn_hat, exponent)
        c_7_k_z = JuMP.@constraint(wm.model, r * lhs_n_z <= dhn / L)

        append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_7_k_y, c_7_k_z])
    end

    for k in 2:pipe_breakpoints-1
        # Add the adjacency constraints for piecewise variables.
        adjacency_p = x_p[a, k-1] + x_p[a, k]
        c_8_k_p = JuMP.@constraint(wm.model, lambda_p[a, k] <= adjacency_p)

        adjacency_n = x_n[a, k-1] + x_n[a, k]
        c_8_k_n = JuMP.@constraint(wm.model, lambda_n[a, k] <= adjacency_n)

        append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_8_k_p, c_8_k_n])
    end
end


########################################## PUMPS ##########################################


"Add constraints associated with modeling a pump's head gain."
function constraint_on_off_pump_head_gain(wm::AbstractPWLRDModel, n::Int, a::Int, node_fr::Int, node_to::Int, q_min_forward::Float64)
    # Get object from the WaterModels reference dictionary.
    pump = ref(wm, n, :pump, a)

    # Gather pump flow, head gain, and status variables.
    qp, g, z = var(wm, n, :qp_pump, a), var(wm, n, :g_pump, a), var(wm, n, :z_pump, a)

    # Gather convex combination variables.
    lambda, x_pw = var(wm, n, :lambda_pump), var(wm, n, :x_pw_pump)

    # Add the required SOS constraints.
    breakpoints = get_pump_flow_lower_breakpoints_positive(pump)
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, k] for k in 1:length(breakpoints)) == z)
    c_2 = JuMP.@constraint(wm.model, sum(x_pw[a, k] for k in 1:length(breakpoints)-1) == z)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, length(breakpoints)] <= x_pw[a, length(breakpoints)-1])

    # Add a constraint for the flow piecewise approximation.
    qp_lhs = sum(breakpoints[k] * lambda[a, k] for k in 1:length(breakpoints))
    c_5 = JuMP.@constraint(wm.model, qp_lhs == qp)

    head_curve_function = _calc_head_curve_function(ref(wm, n, :pump, a))
    head_curve_derivative = _calc_head_curve_derivative(ref(wm, n, :pump, a))

    # Add a constraint that lower-bounds the head gain variable.
    f_all = head_curve_function.(collect(breakpoints))
    gain_lb_expr = sum(f_all[k] .* lambda[a, k] for k in 1:length(breakpoints))
    c_6 = JuMP.@constraint(wm.model, gain_lb_expr <= g)

    # Append the constraint array.
    append!(con(wm, n, :on_off_pump_head_gain, a), [c_1, c_2, c_3, c_4, c_5, c_6])

    for qp_hat in get_pump_flow_upper_breakpoints_positive(pump)
        # Add head gain outer (i.e., upper) approximations.
        f, df = head_curve_function(qp_hat), head_curve_derivative(qp_hat)
        rhs = f * z + df * (qp - qp_hat * z)
        c = JuMP.@constraint(wm.model, g <= rhs)
        append!(con(wm, n, :on_off_pump_head_gain, a), [c])
    end

    for k in 2:length(breakpoints)-1
        # Add the adjacency constraints for piecewise variables.
        adjacency = x_pw[a, k-1] + x_pw[a, k]
        c = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :on_off_pump_head_gain, a), [c])
    end
end


function constraint_on_off_pump_power(wm::AbstractPWLRDModel, n::Int, a::Int, q_min_forward::Float64)
    # Get the number of breakpoints for the pump.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 2)

    # Gather pump flow, power, and status variables.
    q, P = var(wm, n, :qp_pump, a), var(wm, n, :P_pump, a)
    z, lambda = var(wm, n, :z_pump, a), var(wm, n, :lambda_pump)

    # Generate a set of uniform flow breakpoints.
    breakpoints = range(q_min_forward, stop = JuMP.upper_bound(q), length = pump_breakpoints)

    # Add a constraint that lower-bounds the power variable.
    f_ua = _calc_pump_power_ua(wm, n, a, collect(breakpoints))
    power_lb_expr = sum(f_ua[k] * lambda[a, k] for k in 1:pump_breakpoints)
    c_1 = JuMP.@constraint(wm.model, power_lb_expr <= P)

    # Add a constraint that upper-bounds the power variable.
    f_oa = _calc_pump_power_oa(wm, n, a, collect(breakpoints))
    power_ub_expr = sum(f_oa[k] * lambda[a, k] for k in 1:pump_breakpoints)
    c_2 = JuMP.@constraint(wm.model, P <= power_ub_expr)

    # Append the :on_off_pump_power constraint array.
    append!(con(wm, n, :on_off_pump_power)[a], [c_1, c_2])
end