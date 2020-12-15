# Define common PWLRD (piecewise linear relaxation- and direction-based) implementations of water
# distribution network constraints, which use directed flow variables.


########################################## VARIABLES ##########################################


"Create flow-related variables common to all directed flow models for node-connecting components."
function variable_flow(wm::AbstractPWLRDModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    for name in ["pipe", "pump", "regulator", "short_pipe", "valve"]
        # Create directed flow (`qp` and `qn`) variables for each component.
        _variable_component_flow(wm, name; nw=nw, bounded=bounded, report=report)

        # Create directed head difference (`dhp` and `dhn`) variables for each component.
        _variable_component_head_difference(wm, name; nw=nw, bounded=bounded, report=report)

        # Create directed flow binary direction variables (`y`) for each component.
        _variable_component_direction(wm, name; nw=nw, report=report)
    end

    # Create flow-related variables for design components.
    variable_flow_des(wm; nw=nw, bounded=bounded, report=report)

    # Get the number of breakpoints for pumps.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 2)

    # Create weights involved in convex combination constraints for pumps.
    lambda_pump = var(wm, nw)[:lambda_pump] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pump), k in 1:pump_breakpoints],
        base_name="$(nw)_lambda", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :pump, a), "lambda_start", k))

    # Create binary variables involved in convex combination constraints for pumps.
    x_pw_pump = var(wm, nw)[:x_pw_pump] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pump), k in 1:pump_breakpoints-1], base_name="$(nw)_x_pw",
        binary=true, start=comp_start_value(ref(wm, nw, :pump, a), "x_pw_start"))

    # Get the number of breakpoints for pipes.
    pipe_breakpoints = get(wm.ext, :pipe_breakpoints, 2)

    # Create weights involved in convex combination constraints for pipes.
    lambda_p_pipe = var(wm, nw)[:lambda_p_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:pipe_breakpoints],
        base_name="$(nw)_lambda_p", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :pipe, a), "lambda_p_start", k))

    # Create weights involved in convex combination constraints for pipes.
    lambda_n_pipe = var(wm, nw)[:lambda_n_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:pipe_breakpoints],
        base_name="$(nw)_lambda_n", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :pipe, a), "lambda_n_start", k))
    
    # Create binary variables involved in convex combination constraints for pipes.
    x_p_pipe = var(wm, nw)[:x_p_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:pipe_breakpoints-1], base_name="$(nw)_x_p",
        binary=true, start=comp_start_value(ref(wm, nw, :pipe, a), "x_p_start"))

    # Create binary variables involved in convex combination constraints for pipes.
    x_n_pipe = var(wm, nw)[:x_n_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:pipe_breakpoints-1], base_name="$(nw)_x_n",
        binary=true, start=comp_start_value(ref(wm, nw, :pipe, a), "x_n_start"))
end


########################################## PIPES ##########################################


function constraint_pipe_head_loss(wm::PWLRDWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64, L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Get the number of breakpoints for the pipe.
    pipe_breakpoints = get(wm.ext, :pipe_breakpoints, 2)

    # Gather directed pipe flow and direction variables.
    qp, qn, y = var(wm, n, :qp_pipe, a), var(wm, n, :qn_pipe, a), var(wm, n, :y_pipe, a)

    # Gather directed head loss variables.
    dhp, dhn = var(wm, n, :dhp_pipe, a), var(wm, n, :dhn_pipe, a)

    # Gather convex combination variables.
    lambda_p, x_p = var(wm, n, :lambda_p_pipe), var(wm, n, :x_p_pipe)
    lambda_n, x_n = var(wm, n, :lambda_n_pipe), var(wm, n, :x_n_pipe)

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda_p[a, :]) == y)
    c_2 = JuMP.@constraint(wm.model, sum(x_p[a, :]) == y)
    c_3 = JuMP.@constraint(wm.model, lambda_p[a, 1] <= x_p[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda_p[a, end] <= x_p[a, end])

    c_5 = JuMP.@constraint(wm.model, sum(lambda_n[a, :]) == 1.0 - y)
    c_6 = JuMP.@constraint(wm.model, sum(x_n[a, :]) == 1.0 - y)
    c_7 = JuMP.@constraint(wm.model, lambda_n[a, 1] <= x_n[a, 1])
    c_8 = JuMP.@constraint(wm.model, lambda_n[a, end] <= x_n[a, end])

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

    ## Append the constraint array.
    #append!(con(wm, n, :on_off_pipe_head_gain, a), [c_1, c_2, c_3, c_4, c_5, c_6])

    for qp_hat in breakpoints_p
        # Add head loss outer (i.e., lower) approximations.
        lhs_p = _get_head_loss_oa_binary(qp, y, qp_hat, exponent)
        c_7_k = JuMP.@constraint(wm.model, r * lhs_p <= inv(L) * dhp)
        #append!(con(wm, n, :on_off_pipe_head_gain, a), [c_7_k])
    end

    for qn_hat in breakpoints_n
        # Add head loss outer (i.e., lower) approximations.
        lhs_n = _get_head_loss_oa_binary(qn, 1.0 - y, qn_hat, exponent)
        c_7_k = JuMP.@constraint(wm.model, r * lhs_n <= inv(L) * dhn)
        #append!(con(wm, n, :on_off_pipe_head_gain, a), [c_7_k])
    end

    for k in 2:pipe_breakpoints-1
        # Add the adjacency constraints for piecewise variables.
        adjacency_p = x_p[a, k-1] + x_p[a, k]
        c_8_k = JuMP.@constraint(wm.model, lambda_p[a, k] <= adjacency_p)

        adjacency_n = x_n[a, k-1] + x_n[a, k]
        c_8_k = JuMP.@constraint(wm.model, lambda_n[a, k] <= adjacency_n)

        #append!(con(wm, n, :on_off_pipe_head_gain, a), [c_8_k])
    end
end


function constraint_on_off_pipe_head_loss_des(wm::PWLRDWaterModel, n::Int, a::Int, exponent::Float64, node_fr::Int, node_to::Int, L::Float64, resistances)
    # If the number of breakpoints is not positive, no constraints are added.
    num_breakpoints = get(wm.ext, :pipe_breakpoints, 1)

    # Get design pipe direction and status variable references.
    y, z = var(wm, n, :y_des_pipe, a), var(wm, n, :z_des_pipe, a)

    # Get head difference variables for the pipe.
    dhp, dhn = var(wm, n, :dhp_des_pipe, a), var(wm, n, :dhn_des_pipe, a)
    qp, qn = var(wm, n, :qp_des_pipe, a), var(wm, n, :qn_des_pipe, a)

    for (r_id, r) in enumerate(resistances)
        # Get directed flow variables and associated data.
        qp_ub, qn_ub = JuMP.upper_bound(qp[r_id]), JuMP.upper_bound(qn[r_id])

        # Loop over breakpoints strictly between the lower and upper variable bounds.
        for pt in range(0.0, stop = qp_ub, length = num_breakpoints+2)[2:end-1]
            lhs = r * _get_head_loss_oa_binary(qp[r_id], z[r_id], pt, exponent)
            c_1 = JuMP.@constraint(wm.model, lhs <= inv(L) * dhp)
            append!(con(wm, n, :on_off_pipe_head_loss_des)[a], [c_1])
        end

        # Loop over breakpoints strictly between the lower and upper variable bounds.
        for pt in range(0.0, stop = qn_ub, length = num_breakpoints+2)[2:end-1]
            lhs = r * _get_head_loss_oa_binary(qn[r_id], z[r_id], pt, exponent)
            c_2 = JuMP.@constraint(wm.model, lhs <= inv(L) * dhn)
            append!(con(wm, n, :on_off_pipe_head_loss_des)[a], [c_2])
        end
    end

    # Add linear upper bounds for the positive portion of head loss.
    qp_ub = JuMP.upper_bound.(qp)
    slopes_p = resistances .* qp_ub.^(exponent - 1.0)
    c_3 = JuMP.@constraint(wm.model, inv(L)*dhp <= sum(slopes_p .* qp))

    # Add linear upper bounds for the negative portion of head loss.
    qn_ub = JuMP.upper_bound.(qn)
    slopes_n = resistances .* qn_ub.^(exponent - 1.0)
    c_4 = JuMP.@constraint(wm.model, inv(L)*dhn <= sum(slopes_n .* qn))

    # Append the :on_off_pipe_head_loss_des constraint array.
    append!(con(wm, n, :on_off_pipe_head_loss_des)[a], [c_3, c_4])
end


########################################## PUMPS ##########################################


"Add constraints associated with modeling a pump's head gain."
function constraint_on_off_pump_head_gain(wm::PWLRDWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, pc::Array{Float64}, q_min_forward::Float64)
    # Get the number of breakpoints for the pump.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 2)

    # Gather pump flow, head gain, and status variables.
    qp, g, z = var(wm, n, :qp_pump, a), var(wm, n, :g_pump, a), var(wm, n, :z_pump, a)

    # Gather convex combination variables.
    lambda, x_pw = var(wm, n, :lambda_pump), var(wm, n, :x_pw_pump)

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == z)
    c_2 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == z)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

    # Add a constraint for the flow piecewise approximation.
    breakpoints = range(q_min_forward, stop=JuMP.upper_bound(qp), length=pump_breakpoints)
    qp_lhs = sum(breakpoints[k] * lambda[a, k] for k in 1:pump_breakpoints)
    c_5 = JuMP.@constraint(wm.model, qp_lhs == qp)

    head_curve_function = _calc_head_curve_function(ref(wm, n, :pump, a))
    head_curve_derivative = _calc_head_curve_derivative(ref(wm, n, :pump, a))

    # Add a constraint that lower-bounds the head gain variable.
    f_all = head_curve_function.(collect(breakpoints))
    gain_lb_expr = sum(f_all[k] .* lambda[a, k] for k in 1:pump_breakpoints)
    c_6 = JuMP.@constraint(wm.model, gain_lb_expr <= g)

    # Append the constraint array.
    append!(con(wm, n, :on_off_pump_head_gain, a), [c_1, c_2, c_3, c_4, c_5, c_6])

    for qp_hat in breakpoints
        # Add head gain outer (i.e., upper) approximations.
        f, df = head_curve_function(qp_hat), head_curve_derivative(qp_hat)
        rhs = f * z + df * (qp - qp_hat * z)
        c_7_k = JuMP.@constraint(wm.model, g <= rhs)
        append!(con(wm, n, :on_off_pump_head_gain, a), [c_7_k])
    end

    for k in 2:pump_breakpoints-1
        # Add the adjacency constraints for piecewise variables.
        adjacency = x_pw[a, k-1] + x_pw[a, k]
        c_8_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :on_off_pump_head_gain, a), [c_8_k])
    end
end


######################################## OBJECTIVES ########################################


"Instantiate the objective associated with the Optimal Water Flow problem."
function objective_owf_default(wm::PWLRDWaterModel)
    # Get the number of breakpoints for the pump.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 2)

    # Initialize the objective function.
    objective = JuMP.AffExpr(0.0)

    for (n, nw_ref) in nws(wm)
        # Get common variables.
        lambda = var(wm, n, :lambda_pump)

        for (a, pump) in nw_ref[:pump]
            if haskey(pump, "energy_price")
                # Get pump energy price data.
                price = pump["energy_price"]

                # Get pump flow and status variables.
                qp, z = var(wm, n, :qp_pump, a), var(wm, n, :z_pump, a)
                q_min_forward = max(get(pump, "flow_min_forward", _FLOW_MIN), _FLOW_MIN)

                # Generate a set of uniform flow and cubic function breakpoints.
                breakpoints = range(q_min_forward, stop = JuMP.upper_bound(qp), length = pump_breakpoints)
                energy = _calc_pump_energy_ua(wm, n, a, collect(breakpoints))

                # Add the cost corresponding to the current pump's operation.
                cost = sum(price * energy[k] * lambda[a, k] for k in 1:pump_breakpoints)
                JuMP.add_to_expression!(objective, cost)
            else
                Memento.error(_LOGGER, "No cost given for pump \"$(pump["name"])\"")
            end
        end
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end
