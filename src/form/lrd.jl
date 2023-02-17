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
    qp = var(wm, n, :qp_pipe, a)
    dhp = var(wm, n, :dhp_pipe, a)

    # Get the corresponding positive flow partitioning.
    partition_p = get_pipe_flow_partition_positive(ref(wm, n, :pipe, a))

    # Loop over consequential points (i.e., those that have nonzero head loss).
    for flow_value in filter(x -> x > 0.0, partition_p)
        # Add a linear outer approximation of the convex relaxation at `flow_value`.
        lhs = r * _calc_head_loss_oa(qp, y, flow_value, exponent)

        if minimum(abs.(lhs.terms.vals)) >= 1.0e-4
            # Compute a scaling factor to normalize the constraint.
            scalar = _get_scaling_factor(vcat(lhs.terms.vals, [1.0 / L]))

            # Add outer-approximation of the head loss constraint.
            c = JuMP.@constraint(wm.model, scalar * lhs <= scalar * dhp / L)

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

        # Compute a scaling factor to normalize the constraint.
        scalar = _get_scaling_factor(vcat(f_lb_line.terms.vals, [1.0 / L]))

        # Add upper-bounding line of the head loss constraint.
        c = JuMP.@constraint(wm.model, scalar * dhp / L <= scalar * f_lb_line)

        # Append the :pipe_head_loss constraint array.
        append!(con(wm, n, :pipe_head_loss)[a], [c])
    elseif qp_max == 0.0
        c = JuMP.@constraint(wm.model, dhp == 0.0)
        append!(con(wm, n, :pipe_head_loss)[a], [c])
    else
        f_q = r * qp_max^exponent
        scalar = _get_scaling_factor([f_q == 0.0 ? 1.0 : f_q, 1.0 / L])
        c = JuMP.@constraint(wm.model, scalar * dhp / L == scalar * f_q * y)
        append!(con(wm, n, :pipe_head_loss)[a], [c])
    end

    # Get variables for negative flow and head difference.
    qn = var(wm, n, :qn_pipe, a)
    dhn = var(wm, n, :dhn_pipe, a)

    # Get the corresponding negative flow partitioning (negated).
    partition_n = sort(-get_pipe_flow_partition_negative(ref(wm, n, :pipe, a)))

    # Loop over consequential points (i.e., those that have nonzero head loss).
    for flow_value in filter(x -> x > 0.0, partition_n)
        # Add a linear outer approximation of the convex relaxation at `flow_value`.
        lhs = r * _calc_head_loss_oa(qn, 1.0 - y, flow_value, exponent)

        if minimum(abs.(lhs.terms.vals)) >= 1.0e-4
            # Compute a scaling factor to normalize the constraint.
            scalar = _get_scaling_factor(vcat(lhs.terms.vals, [1.0 / L]))

            # Add outer-approximation of the head loss constraint.
            c = JuMP.@constraint(wm.model, scalar * lhs <= scalar * dhn / L)

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

        # Compute a scaling factor to normalize the constraint.
        scalar = _get_scaling_factor(vcat(f_lb_line.terms.vals, [1.0 / L]))

        # Add upper-bounding line of the head loss constraint.
        c = JuMP.@constraint(wm.model, scalar * dhn / L <= scalar * f_lb_line)

        # Append the :pipe_head_loss constraint array.
        append!(con(wm, n, :pipe_head_loss)[a], [c])
    elseif qn_max == 0.0
        c = JuMP.@constraint(wm.model, dhn == 0.0)
        append!(con(wm, n, :pipe_head_loss)[a], [c])
    else
        f_q = r * qn_max^exponent
        scalar = _get_scaling_factor([f_q == 0.0 ? 1.0 : f_q, 1.0 / L])
        c = JuMP.@constraint(wm.model, scalar * dhn / L == scalar * f_q * (1.0 - y))
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
        scalar = _get_scaling_factor(vcat(lhs_1.terms.vals, [1.0 / L]))
        c_1 = JuMP.@constraint(wm.model, scalar * lhs_1 <= scalar * dhp / L)

        scalar = _get_scaling_factor(vcat(lhs_2.terms.vals, [1.0 / L]))
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
        scalar = _get_scaling_factor(vcat(f_lb_line_z.terms.vals, [1.0 / L]))
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
        scalar = _get_scaling_factor(vcat(lhs_1.terms.vals, [1.0 / L]))
        c_1 = JuMP.@constraint(wm.model, scalar * lhs_1 <= scalar * dhn / L)

        scalar = _get_scaling_factor(vcat(lhs_2.terms.vals, [1.0 / L]))
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
        scalar = _get_scaling_factor(vcat(f_lb_line_z.terms.vals, [1.0 / L]))
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
    if(n==1)
        println("Running LRD pump head gain")
    end
    # Get variables for positive flow, head difference, and pump status.
    qp = var(wm, n, :qp_pump, a)
    g = var(wm, n, :g_pump, a)
    z = var(wm, n, :z_pump, a)

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
            # Compute a scaling factor to normalize the constraint.
            scalar = _get_scaling_factor([1.0, f, df * flow_value])

            # Add the outer-approximation constraint for the pump.
            rhs = f * z + df * (qp - flow_value * z)
            c_1 = JuMP.@constraint(wm.model, scalar * g <= scalar * rhs)

            # Append the :on_off_pump_head_gain constraint array.
            append!(con(wm, n, :on_off_pump_head_gain)[a], [c_1])
        end
    end

    if qp_min < qp_max
        # Collect relevant expressions for the lower-bounding line.
        g_1, g_2 = head_curve_func(qp_min), head_curve_func(qp_max)
        g_slope = (g_2 - g_1) / (qp_max - qp_min)
        g_lb_line = g_slope * (qp - qp_min * z) + g_1 * z

        # Compute a scaling factor to normalize the constraint.
        scalar = _get_scaling_factor(vcat([1.0], g_lb_line.terms.vals))

        # Add the lower-bounding line for the head gain curve.
        c_2 = JuMP.@constraint(wm.model, scalar * g_lb_line <= scalar * g)

        # Append the :on_off_pump_head_gain constraint array.
        append!(con(wm, n, :on_off_pump_head_gain)[a], [c_2])
    end
end


"Add constraints associated with modeling an expansion pump's head gain."
function constraint_on_off_pump_head_gain_ne(
    wm::AbstractLRDModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    q_min_forward::Float64,
)
    # Get variables for positive flow, head difference, and pump status.
    qp = var(wm, n, :qp_ne_pump, a)
    g = var(wm, n, :g_ne_pump, a)
    z = var(wm, n, :z_ne_pump, a)

    # Calculate the head curve function and its derivative.
    # Get pump head curve function and its derivative.
    head_curve_func = ref(wm, n, :ne_pump, a, "head_curve_function")
    head_curve_deriv = ref(wm, n, :ne_pump, a, "head_curve_derivative")
    partition = get_pump_flow_partition(ref(wm, n, :ne_pump, a))
    qp_min, qp_max = minimum(partition), maximum(partition)

    # Loop over partition points strictly between the lower and upper variable bounds.
    for flow_value in partition
        # Compute head gain and derivative at the point.
        f, df = head_curve_func(flow_value), head_curve_deriv(flow_value)

        if abs(df) >= 1.0e-4 # Only add an outer-approximation if the derivative isn't too small.
            # Compute a scaling factor to normalize the constraint.
            scalar = _get_scaling_factor([1.0, f, df * flow_value])

            # Add the outer-approximation constraint for the pump.
            rhs = f * z + df * (qp - flow_value * z)
            c_1 = JuMP.@constraint(wm.model, scalar * g <= scalar * rhs)

            # Append the :on_off_pump_head_gain constraint array.
            append!(con(wm, n, :on_off_pump_head_gain_ne)[a], [c_1])
        end
    end

    if qp_min < qp_max
        # Collect relevant expressions for the lower-bounding line.
        g_1, g_2 = head_curve_func(qp_min), head_curve_func(qp_max)
        g_slope = (g_2 - g_1) / (qp_max - qp_min)
        g_lb_line = g_slope * (qp - qp_min * z) + g_1 * z

        # Compute a scaling factor to normalize the constraint.
        scalar = _get_scaling_factor(vcat([1.0], g_lb_line.terms.vals))

        # Add the lower-bounding line for the head gain curve.
        c_2 = JuMP.@constraint(wm.model, scalar * g_lb_line <= scalar * g)

        # Append the :on_off_pump_head_gain constraint array.
        append!(con(wm, n, :on_off_pump_head_gain_ne)[a], [c_2])
    end
end
