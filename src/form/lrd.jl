# Define common LRD (linear relaxation- and direction-based) implementations of water
# distribution network constraints, which use directed flow variables.


########################################## PIPES ##########################################

function constraint_pipe_head_loss(
    wm::AbstractLRDModel,
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
    # Get the variable for flow directionality.
    y = var(wm, n, :y_pipe, a)

    # Get variables for positive flow and head difference.
    qp, dhp = var(wm, n, :qp_pipe, a), var(wm, n, :dhp_pipe, a)

    # Get the corresponding positive flow partitioning.
    partition_p = get_pipe_flow_partition_positive(ref(wm, n, :pipe, a))

    # Loop over consequential points (i.e., those that have nonzero head loss).
    for flow_value in filter(x -> x > 0.0, partition_p)
        # Add a linear outer approximation of the convex relaxation at `flow_value`.
        lhs = r * _calc_head_loss_oa(qp, y, flow_value, exponent)

        if minimum(abs.(lhs.terms.vals)) >= 1.0e-4
            # Add outer-approximation of the head loss constraint.
            c = JuMP.@constraint(wm.model, lhs <= dhp / L)

            # Append the :pipe_head_loss constraint array.
            append!(con(wm, n, :pipe_head_loss)[a], [c])
        end
    end

    # Get the corresponding min/max positive directed flows (when active).
    qp_min_forward, qp_max = max(0.0, q_min_forward), maximum(partition_p)

    if qp_min_forward != qp_max
        # Compute scaled version of the head loss overapproximation.
        f_1, f_2 = r * qp_min_forward^exponent, r * qp_max^exponent
        f_slope = (f_2 - f_1) / (qp_max - qp_min_forward)
        f_lb_line = f_slope * (qp - qp_min_forward * y) + f_1 * y

        # Add upper-bounding line of the head loss constraint.
        c = JuMP.@constraint(wm.model, dhp / L <= f_lb_line)

        # Append the :pipe_head_loss constraint array.
        append!(con(wm, n, :pipe_head_loss)[a], [c])
    end

    # Get variables for negative flow and head difference.
    qn, dhn = var(wm, n, :qn_pipe, a), var(wm, n, :dhn_pipe, a)

    # Get the corresponding negative flow partitioning (negated).
    partition_n = sort(-get_pipe_flow_partition_negative(ref(wm, n, :pipe, a)))

    # Loop over consequential points (i.e., those that have nonzero head loss).
    for flow_value in filter(x -> x > 0.0, partition_n)
        # Add a linear outer approximation of the convex relaxation at `flow_value`.
        lhs = r * _calc_head_loss_oa(qn, 1.0 - y, flow_value, exponent)

        if minimum(abs.(lhs.terms.vals)) >= 1.0e-4
            # Add outer-approximation of the head loss constraint.
            c = JuMP.@constraint(wm.model, lhs <= dhn / L)

            # Append the :pipe_head_loss constraint array.
            append!(con(wm, n, :pipe_head_loss)[a], [c])
        end
    end

    # Get the corresponding maximum negative directed flow (when active).
    qn_min_forward, qn_max = max(0.0, -q_max_reverse), maximum(partition_n)

    if qn_min_forward != qn_max
        # Compute scaled version of the head loss overapproximation.
        f_1, f_2 = r * qn_min_forward^exponent, r * qn_max^exponent
        f_slope = (f_2 - f_1) / (qn_max - qn_min_forward)
        f_lb_line = f_slope * (qn - qn_min_forward * (1.0 - y)) + f_1 * (1.0 - y)

        # Add upper-bounding line of the head loss constraint.
        c = JuMP.@constraint(wm.model, dhn / L <= f_lb_line)

        # Append the :pipe_head_loss constraint array.
        append!(con(wm, n, :pipe_head_loss)[a], [c])
    end
end


function constraint_on_off_des_pipe_head_loss(
    wm::AbstractLRDModel,
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
    # Get the variable for flow directionality.
    y = var(wm, n, :y_des_pipe, a)
    z = var(wm, n, :z_des_pipe, a)

    # Get variables for positive flow and head difference.
    qp = var(wm, n, :qp_des_pipe, a)
    dhp = var(wm, n, :dhp_des_pipe, a)

    # Get the corresponding positive flow partitioning.
    des_pipe = ref(wm, n, :des_pipe, a)
    partition_p = get_pipe_flow_partition_positive(des_pipe)

    # Loop over consequential points (i.e., those that have nonzero head loss).
    for flow_value in filter(x -> x > 0.0, partition_p)
        # Get linear outer approximations of the convex relaxation at `flow_value`.
        lhs_1 = r * _calc_head_loss_oa(qp, z, flow_value, exponent)
        lhs_2 = r * _calc_head_loss_oa(qp, y, flow_value, exponent)

        # Add outer approximations of the head loss constraint.
        min_val = minimum(abs.(filter(x -> x != 0.0, lhs_1.terms.vals)))
        scalar = 10^(-0.25 * log10(min_val)) # Helps with scaling issues.

        c_1 = JuMP.@constraint(wm.model, scalar * lhs_1 <= scalar * dhp / L)
        c_2 = JuMP.@constraint(wm.model, scalar * lhs_2 <= scalar * dhp / L)

        # Append the :on_off_des_pipe_head_loss constraint array.
        append!(con(wm, n, :on_off_des_pipe_head_loss)[a], [c_1, c_2])
    end

    # Get the corresponding min/max positive directed flows (when active).
    qp_min_forward, qp_max = max(0.0, q_min_forward), maximum(partition_p)

    if qp_min_forward != qp_max
        # Compute scaled versions of the head loss overapproximations.
        f_1, f_2 = r * qp_min_forward^exponent, r * qp_max^exponent
        f_slope = (f_2 - f_1) / (qp_max - qp_min_forward)
        f_lb_line_z = f_slope * (qp - qp_min_forward * z) + f_1 * z

        # Add upper-bounding lines of the head loss constraint.
        min_val = minimum(abs.(filter(x -> x != 0.0, f_lb_line_z.terms.vals)))
        scalar = 10^(-0.25 * log10(min_val)) # Helps with scaling issues.
        c = JuMP.@constraint(wm.model, scalar * dhp / L <= scalar * f_lb_line_z)

        # Append the :on_off_des_pipe_head_loss constraint array.
        append!(con(wm, n, :on_off_des_pipe_head_loss)[a], [c])
    end

    # Get variables for negative flow and head difference.
    qn, dhn = var(wm, n, :qn_des_pipe, a), var(wm, n, :dhn_des_pipe, a)

    # Get the corresponding negative flow partitioning (negated).
    partition_n = sort(-get_pipe_flow_partition_negative(ref(wm, n, :des_pipe, a)))

    # Loop over consequential points (i.e., those that have nonzero head loss).
    for flow_value in filter(x -> x > 0.0, partition_n)
        # Get linear outer approximations of the convex relaxation at `flow_value`.
        lhs_1 = r * _calc_head_loss_oa(qn, z, flow_value, exponent)
        lhs_2 = r * _calc_head_loss_oa(qn, 1.0 - y, flow_value, exponent)

        # Add outer approximations of the head loss constraint.
        min_val = minimum(abs.(filter(x -> x != 0.0, lhs_1.terms.vals)))
        scalar = 10^(-0.25 * log10(min_val)) # Helps with scaling issues.
        c_1 = JuMP.@constraint(wm.model, scalar * lhs_1 <= scalar * dhn / L)
        c_2 = JuMP.@constraint(wm.model, scalar * lhs_2 <= scalar * dhn / L)

        # Append the :on_off_des_pipe_head_loss constraint array.
        append!(con(wm, n, :on_off_des_pipe_head_loss)[a], [c_1, c_2])
    end

    # Get the corresponding maximum negative directed flow (when active).
    qn_min_forward, qn_max = max(0.0, -q_max_reverse), maximum(partition_n)

    if qn_min_forward != qn_max
        # Compute scaled versions of the head loss overapproximations.
        f_1, f_2 = r * qn_min_forward^exponent, r * qn_max^exponent
        f_slope = (f_2 - f_1) / (qn_max - qn_min_forward)
        f_lb_line_z = f_slope * (qn - qn_min_forward * z) + f_1 * z

        # Add upper-bounding lines of the head loss constraint.
        min_val = minimum(abs.(filter(x -> x != 0.0, f_lb_line_z.terms.vals)))
        scalar = 10^(-0.25 * log10(min_val)) # Helps with scaling issues.
        c = JuMP.@constraint(wm.model, scalar * dhn / L <= scalar * f_lb_line_z)

        # Append the :on_off_des_pipe_head_loss constraint array.
        append!(con(wm, n, :on_off_des_pipe_head_loss)[a], [c])
    end
end


########################################## PUMPS ##########################################

"Add constraints associated with modeling a pump's head gain."
function constraint_on_off_pump_head_gain(
    wm::AbstractLRDModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    q_min_forward::Float64,
)
    # Get variables for positive flow, head difference, and pump status.
    qp, g, z = var(wm, n, :qp_pump, a), var(wm, n, :g_pump, a), var(wm, n, :z_pump, a)

    # Calculate the head curve function and its derivative.
    # Get pump head curve function and its derivative.
    head_curve_func = ref(wm, n, :pump, a, "head_curve_function")
    head_curve_deriv = ref(wm, n, :pump, a, "head_curve_derivative")
    partition = get_pump_flow_partition(ref(wm, n, :pump, a))
    qp_min, qp_max = minimum(partition), maximum(partition)

    # Loop over partition points strictly between the lower and upper variable bounds.
    for flow_value in partition
        # Compute head gain and derivative at the point.
        f, df = head_curve_func(flow_value), head_curve_deriv(flow_value)

        if abs(df) >= 1.0e-4 # Only add an outer-approximation if the derivative isn't too small.
            # Add the outer-approximation constraint for the pump.
            c_1 = JuMP.@constraint(wm.model, g <= f * z + df * (qp - flow_value * z))
    
            # Append the :on_off_pump_head_gain constraint array.
            append!(con(wm, n, :on_off_pump_head_gain)[a], [c_1])
        end
    end

    if qp_min < qp_max
        # Collect relevant expressions for the lower-bounding line.
        g_1, g_2 = head_curve_func(qp_min), head_curve_func(qp_max)
        g_slope = (g_2 - g_1) / (qp_max - qp_min)
        g_lb_line = g_slope * (qp - qp_min * z) + g_1 * z

        # Add the lower-bounding line for the head gain curve.
        c_2 = JuMP.@constraint(wm.model, g_lb_line <= g)

        # Append the :on_off_pump_head_gain constraint array.
        append!(con(wm, n, :on_off_pump_head_gain)[a], [c_2])
    end
end
