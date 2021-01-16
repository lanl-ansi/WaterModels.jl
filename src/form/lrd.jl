# Define common LRD (linear relaxation- and direction-based) implementations of water
# distribution network constraints, which use directed flow variables.


########################################## PIPES ##########################################

function constraint_pipe_head_loss(
    wm::AbstractLRDModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
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
        lhs = _calc_head_loss_oa(qp, y, pt, exponent)

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
        lhs = _calc_head_loss_oa(qn, 1.0 - y, pt, exponent)

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
    y, z = var(wm, n, :y_des_pipe, a), var(wm, n, :z_des_pipe, a)

    # Get variables for positive flow and head difference.
    qp, dhp = var(wm, n, :qp_des_pipe, a), var(wm, n, :dhp_des_pipe, a)
    qp_min_forward, qp_ub = max(0.0, q_min_forward), JuMP.upper_bound(qp)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(qp_min_forward, stop = JuMP.upper_bound(qp), length = num_breakpoints+2)[2:end-1]
        # Get a linear outer approximation of the convex relaxation at `pt`.
        lhs = _calc_head_loss_oa(qp, z, pt, exponent)

        # Add outer-approximation of the head loss constraint.
        c = JuMP.@constraint(wm.model, r * lhs <= dhp / L)

        # Append the :on_off_des_pipe_head_loss constraint array.
        append!(con(wm, n, :on_off_des_pipe_head_loss)[a], [c])
    end

    # Get the upper-bounding line for the qp curve.
    if qp_min_forward < qp_ub
        dhp_1, dhp_2 = r * qp_min_forward^exponent, r * qp_ub^exponent
        dhp_slope = (dhp_2 - dhp_1) * inv(qp_ub - qp_min_forward)
        dhp_ub_line_y = dhp_slope * (qp - qp_min_forward * y) + dhp_1 * y
        dhp_ub_line_z = dhp_slope * (qp - qp_min_forward * z) + dhp_1 * z

        # Add upper-bounding lines of the head loss constraint.
        c_1 = JuMP.@constraint(wm.model, dhp / L <= dhp_ub_line_y)
        c_2 = JuMP.@constraint(wm.model, dhp / L <= dhp_ub_line_z)

        # Append the :on_off_des_pipe_head_loss constraint array.
        append!(con(wm, n, :on_off_des_pipe_head_loss)[a], [c_1, c_2])
    end

    # Get variables for negative flow and head difference.
    qn, dhn = var(wm, n, :qn_des_pipe, a), var(wm, n, :dhn_des_pipe, a)
    qn_min_forward, qn_ub = max(0.0, -q_max_reverse), JuMP.upper_bound(qn)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(qn_min_forward, stop = JuMP.upper_bound(qn), length = num_breakpoints+2)[2:end-1]
        # Get a linear outer approximation of the convex relaxation at `pt`.
        lhs = _calc_head_loss_oa(qn, z, pt, exponent)

        # Add lower-bounding line of the head loss constraint.
        c = JuMP.@constraint(wm.model, r * lhs <= dhn / L)

        # Append the :on_off_des_pipe_head_loss constraint array.
        append!(con(wm, n, :on_off_des_pipe_head_loss)[a], [c])
    end

    # Get the lower-bounding line for the qn curve.
    if qn_min_forward < qn_ub
        dhn_1, dhn_2 = r * qn_min_forward^exponent, r * qn_ub^exponent
        dhn_slope = (dhn_2 - dhn_1) * inv(qn_ub - qn_min_forward)
        dhn_ub_line_y = dhn_slope * (qn - qn_min_forward * (1.0 - y)) + dhn_1 * (1.0 - y)
        dhn_ub_line_z = dhn_slope * (qn - qn_min_forward * z) + dhn_1 * z

        # Add upper-bounding line of the head loss constraint.
        c_1 = JuMP.@constraint(wm.model, dhn / L <= dhn_ub_line_y)
        c_2 = JuMP.@constraint(wm.model, dhn / L <= dhn_ub_line_z)

        # Append the :on_off_des_pipe_head_loss constraint array.
        append!(con(wm, n, :on_off_des_pipe_head_loss)[a], [c_1, c_2])
    end
end


########################################## PUMPS ##########################################

"Add constraints associated with modeling a pump's head gain."
function constraint_on_off_pump_head_gain(wm::AbstractLRDModel, n::Int, a::Int, node_fr::Int, node_to::Int, q_min_forward::Float64)
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