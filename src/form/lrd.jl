# Define common LRD (linear relaxation- and direction-based) implementations of water
# distribution network constraints, which use directed flow variables.


########################################## PIPES ##########################################

function constraint_pipe_head_loss(
    wm::LRDWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Get the number of breakpoints for the pipe.
    num_breakpoints = get(wm.ext, :pipe_breakpoints, 1)

    # Get the variable for flow directionality.
    y = var(wm, n, :y_pipe, a)

    # Get variables for positive flow and head difference.
    qp, dhp = var(wm, n, :qp_pipe, a), var(wm, n, :dhp_pipe, a)
    qp_min_forward, qp_ub = max(0.0, q_min_forward), JuMP.upper_bound(qp)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(qp_min_forward, stop = JuMP.upper_bound(qp), length = num_breakpoints+2)[2:end-1]
        # Add a linear outer approximation of the convex relaxation at `pt`.
        lhs = _get_head_loss_oa_binary(qp, y, pt, exponent)

        # Add outer-approximation of the head loss constraint.
        c = JuMP.@constraint(wm.model, r * lhs - inv(L) * dhp <= 0.0)

        # Append the :pipe_head_loss constraint array.
        append!(con(wm, n, :pipe_head_loss)[a], [c])
    end

    if qp_min_forward < qp_ub
        dhp_1, dhp_2 = r * qp_min_forward^exponent, r * qp_ub^exponent
        dhp_slope = (dhp_2 - dhp_1) * inv(qp_ub - qp_min_forward)
        dhp_lb_line = dhp_slope * (qp - qp_min_forward * y) + dhp_1 * y

        # Add upper-bounding line of the head loss constraint.
        c = JuMP.@constraint(wm.model, inv(L) * dhp - dhp_lb_line <= 0.0)

        # Append the :pipe_head_loss constraint array.
        append!(con(wm, n, :pipe_head_loss)[a], [c])
    end

    # Get variables for negative flow and head difference.
    qn, dhn = var(wm, n, :qn_pipe, a), var(wm, n, :dhn_pipe, a)
    qn_min_forward, qn_ub = max(0.0, -q_max_reverse), JuMP.upper_bound(qn)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(qn_min_forward, stop = JuMP.upper_bound(qn), length = num_breakpoints+2)[2:end-1]
        # Add a linear outer approximation of the convex relaxation at `pt`.
        lhs = _get_head_loss_oa_binary(qn, 1.0 - y, pt, exponent)

        # Add outer-approximation of the head loss constraint.
        c = JuMP.@constraint(wm.model, r * lhs - inv(L) * dhn <= 0.0)

        # Append the :pipe_head_loss constraint array.
        append!(con(wm, n, :pipe_head_loss)[a], [c])
    end

    if qn_min_forward < qn_ub
        dhn_1, dhn_2 = r * qn_min_forward^exponent, r * qn_ub^exponent
        dhn_slope = (dhn_2 - dhn_1) * inv(qn_ub - qn_min_forward)
        dhn_lb_line = dhn_slope * (qn - qn_min_forward * (1.0 - y)) + dhn_1 * (1.0 - y)

        # Add upper-bounding line of the head loss constraint.
        c = JuMP.@constraint(wm.model, inv(L) * dhn - dhn_lb_line <= 0.0)

        # Append the :pipe_head_loss constraint array.
        append!(con(wm, n, :pipe_head_loss)[a], [c])
    end
end


function constraint_on_off_des_pipe_head_loss(
    wm::AbstractLRDModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Get the number of breakpoints for the pipe.
    num_breakpoints = get(wm.ext, :pipe_breakpoints, 1)

    # Get the variable for flow directionality.
    y = var(wm, n, :y_des_pipe, a)

    # Get variables for positive flow and head difference.
    qp, dhp = var(wm, n, :qp_des_pipe, a), var(wm, n, :dhp_des_pipe, a)
    qp_min_forward, qp_ub = max(0.0, q_min_forward), JuMP.upper_bound(qp)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(qp_min_forward, stop = JuMP.upper_bound(qp), length = num_breakpoints+2)[2:end-1]
        # Add a linear outer approximation of the convex relaxation at `pt`.
        lhs = _get_head_loss_oa_binary(qp, y, pt, exponent)

        # Add the normalized constraint to the model.
        c = JuMP.@constraint(wm.model, r * lhs - inv(L) * dhp <= 0.0)

        # Append the :on_off_des_pipe_head_loss constraint array.
        append!(con(wm, n, :on_off_des_pipe_head_loss)[a], [c])
    end

    # Get the lower-bounding line for the qp curve.
    if qp_min_forward < qp_ub
        dhp_1, dhp_2 = r * qp_min_forward^exponent, r * qp_ub^exponent
        dhp_slope = (dhp_2 - dhp_1) * inv(qp_ub - qp_min_forward)
        dhp_lb_line = dhp_slope * (qp - qp_min_forward * y) + dhp_1 * y

        # Add the normalized constraint to the model.
        c = JuMP.@constraint(wm.model, inv(L) * dhp - dhp_lb_line <= 0.0)

        # Append the :on_off_des_pipe_head_loss constraint array.
        append!(con(wm, n, :on_off_des_pipe_head_loss)[a], [c])
    end

    # Get variables for negative flow and head difference.
    qn, dhn = var(wm, n, :qn_des_pipe, a), var(wm, n, :dhn_des_pipe, a)
    qn_min_forward, qn_ub = max(0.0, -q_max_reverse), JuMP.upper_bound(qn)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(qn_min_forward, stop = JuMP.upper_bound(qn), length = num_breakpoints+2)[2:end-1]
        # Add a linear outer approximation of the convex relaxation at `pt`.
        lhs = _get_head_loss_oa_binary(qn, 1.0 - y, pt, exponent)

        # Add the normalized constraint to the model.
        c = JuMP.@constraint(wm.model, r * lhs - inv(L) * dhn <= 0.0)

        # Append the :on_off_des_pipe_head_loss constraint array.
        append!(con(wm, n, :on_off_des_pipe_head_loss)[a], [c])
    end

    # Get the lower-bounding line for the qn curve.
    if qn_min_forward < qn_ub
        dhn_1, dhn_2 = r * qn_min_forward^exponent, r * qn_ub^exponent
        dhn_slope = (dhn_2 - dhn_1) * inv(qn_ub - qn_min_forward)
        dhn_lb_line = dhn_slope * (qn - qn_min_forward * (1.0 - y)) + dhn_1 * (1.0 - y)

        # Add the normalized constraint to the model.
        c = JuMP.@constraint(wm.model, inv(L) * dhn - dhn_lb_line <= 0.0)

        # Append the :on_off_des_pipe_head_loss constraint array.
        append!(con(wm, n, :on_off_des_pipe_head_loss)[a], [c])
    end
end


########################################## PUMPS ##########################################

"Add constraints associated with modeling a pump's head gain."
function constraint_on_off_pump_head_gain(wm::LRDWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, pc::Array{Float64}, q_min_forward::Float64)
    # Get the number of breakpoints for the pump.
    num_breakpoints = get(wm.ext, :pump_breakpoints, 1)

    # Get variables for positive flow, head difference, and pump status.
    qp, g, z = var(wm, n, :qp_pump, a), var(wm, n, :g_pump, a), var(wm, n, :z_pump, a)
    qp_min_forward, qp_ub = max(0.0, q_min_forward), JuMP.upper_bound(qp)

    # Calculate the head curve function and its derivative.
    head_curve_func = _calc_head_curve_function(ref(wm, n, :pump, a))
    head_curve_deriv = _calc_head_curve_derivative(ref(wm, n, :pump, a))

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(qp_min_forward, stop = JuMP.upper_bound(qp), length = num_breakpoints+2)[2:end-1]
        # Compute head gain and derivative at the point.
        f, df = head_curve_func(pt), head_curve_deriv(pt)

        # Add the outer-approximation constraint for the pump.
        c = JuMP.@constraint(wm.model, g <= f * z + df * (qp - pt * z))

        # Append the :on_off_pump_head_gain constraint array.
        append!(con(wm, n, :on_off_pump_head_gain)[a], [c])
    end

    if qp_min_forward < qp_ub
        # Collect relevant expressions for the lower-bounding line.
        g_1, g_2 = head_curve_func(qp_min_forward), head_curve_func(qp_ub)
        g_slope = (g_2 - g_1) * inv(qp_ub - qp_min_forward)
        g_lb_line = g_slope * (qp - qp_min_forward * z) + g_1 * z

        # Add the lower-bounding line for the head gain curve.
        c = JuMP.@constraint(wm.model, g_lb_line <= g)

        # Append the :on_off_pump_head_gain constraint array.
        append!(con(wm, n, :on_off_pump_head_gain)[a], [c])
    end
end


######################################## OBJECTIVES ########################################


"Instantiate the objective associated with the Optimal Water Flow problem."
function objective_owf_default(wm::LRDWaterModel)
    # Initialize the objective function.
    objective = JuMP.AffExpr(0.0)

    for (n, nw_ref) in nws(wm)
        for (a, pump) in nw_ref[:pump]
            # Ensure that the pump has an associated energy price.
            @assert haskey(pump, "energy_price")

            # Get flow-related variables and data.
            qp, z = var(wm, n, :qp_pump, a), var(wm, n, :z_pump, a)
            qp_lb, qp_ub = pump["flow_min_forward"], JuMP.upper_bound(qp)
            f_ua = _calc_pump_energy_ua(wm, n, a, [qp_lb, qp_ub])

            # Build a linear under-approximation of the cost.
            slope = (f_ua[2] - f_ua[1]) * inv(qp_ub - qp_lb)
            energy = (slope * (qp - qp_lb * z) + f_ua[1] * z)
            JuMP.add_to_expression!(objective, pump["energy_price"] * energy)
        end
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end
