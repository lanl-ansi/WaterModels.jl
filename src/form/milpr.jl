# Define MILP-R (relaxation-based mixed-integer linear programming)
# implementations of common water distribution model specifications.

function _get_owf_oa(q::JuMP.VariableRef, z::JuMP.VariableRef, q_hat::Float64, coeffs::Array{Float64})
    f = coeffs[1]*q_hat^3 + coeffs[2]*q_hat^2 + coeffs[3]*q_hat
    df = 3.0*coeffs[1]*q_hat^2 + 2.0*coeffs[2]*q_hat + coeffs[3]
    return f*z + df*(q - q_hat*z)
end


function _get_head_loss_oa(q::JuMP.VariableRef, q_hat::Float64, alpha::Float64)
    return q_hat^alpha + alpha * q_hat^(alpha - 1.0) * (q - q_hat)
end


function _get_head_loss_cv_oa(q::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, q_hat::Float64, alpha::Float64)
    return q_hat^alpha*z + alpha * q_hat^(alpha - 1.0) * (q - q_hat*z)
end


function _get_head_gain_oa(q::JuMP.VariableRef, z::JuMP.VariableRef, q_hat::Float64, curve_fun::Array{Float64})
    f = curve_fun[1]*q_hat^2 + curve_fun[2]*q_hat + curve_fun[3]
    df = 2.0 * curve_fun[1] * q_hat + curve_fun[2]
    return f * z + df * (q - q_hat * z)
end


function variable_pump_operation(wm::AbstractMILPRModel; nw::Int=wm.cnw, report::Bool=true)
    # If the number of breakpoints is not positive, return.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 0)

    if pump_breakpoints > 0
        # Create weights involved in convex combination constraints.
        lambda = var(wm, nw)[:lambda_pump] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pump), k in 1:pump_breakpoints],
            base_name="$(nw)_lambda", lower_bound=0.0, upper_bound=1.0,
            start=comp_start_value(ref(wm, nw, :pump, a), "lambda_start", k))

        # Create binary variables involved in convex combination constraints.
        x_pw = var(wm, nw)[:x_pw_pump] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pump), k in 1:pump_breakpoints-1],
            base_name="$(nw)_x_pw", binary=true,
            start=comp_start_value(ref(wm, nw, :pump, a), "x_pw_start", k))
    end
end


"Pump head gain constraint when the pump status is ambiguous."
function constraint_pump_head_gain(wm::AbstractMILPRModel, n::Int, a::Int, node_fr::Int, node_to::Int, pc::Array{Float64})
    # If the number of breakpoints is not positive, no constraints are added.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 0)
    if pump_breakpoints <= 0 return end

    # Gather flow, head gain, and convex combination variables.
    qp, g, z = var(wm, n, :qp_pump, a), var(wm, n, :g, a), var(wm, n, :z_pump, a)
    lambda, x_pw = [var(wm, n, :lambda_pump), var(wm, n, :x_pw_pump)]

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == z)
    c_2 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == z)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

    # Add a constraint for the flow piecewise approximation.
    breakpoints = range(0.0, stop=JuMP.upper_bound(qp), length=pump_breakpoints)
    qp_lhs = sum(breakpoints[k] * lambda[a, k] for k in 1:pump_breakpoints)
    c_5 = JuMP.@constraint(wm.model, qp_lhs == qp)

    # Append the constraint array.
    append!(con(wm, n, :head_gain, a), [c_1, c_2, c_3, c_4, c_5])

    for k in 2:pump_breakpoints-1
        # Add the adjacency constraints for piecewise variables.
        adjacency = x_pw[a, k-1] + x_pw[a, k]
        c_6_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :head_gain, a), [c_6_k])
    end

    # Add head gain outer (i.e., upper) approximations.
    for qp_hat in breakpoints
        lhs = _get_head_gain_oa(qp, z, qp_hat, pc)
        c_7_k = JuMP.@constraint(wm.model, g <= lhs)
        append!(con(wm, n, :head_gain, a), [c_7_k])
    end
end


function constraint_pipe_head_loss_des(wm::AbstractMILPRModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, resistances)
    # If the number of breakpoints is not positive, no constraints are added.
    pipe_breakpoints = get(wm.ext, :pipe_breakpoints, 0)
    if pipe_breakpoints <= 0 return end

    for (r_id, r) in enumerate(resistances)
        # Add constraints corresponding to positive outer approximations.
        qp, dhp = var(wm, n, :qp_des, a)[r_id], var(wm, n, :dhp_pipe, a)
        qp_ub = JuMP.upper_bound(qp)

        for qp_hat in range(0.0, stop=qp_ub, length=pipe_breakpoints)
            lhs = r * _get_head_loss_oa(qp, qp_hat, alpha)
            c_p = JuMP.@constraint(wm.model, lhs <= inv(L) * dhp)
            append!(con(wm, n, :head_loss)[a], [c_p])
        end

        # Add constraints corresponding to negative outer approximations.
        qn, dhn = var(wm, n, :qn_des, a)[r_id], var(wm, n, :dhn_pipe, a)
        qn_ub = JuMP.upper_bound(qn)

        for qn_hat in range(0.0, stop=qn_ub, length=pipe_breakpoints)
            lhs = r * _get_head_loss_oa(qn, qn_hat, alpha)
            c_n = JuMP.@constraint(wm.model, lhs <= inv(L) * dhn)
            append!(con(wm, n, :head_loss)[a], [c_n])
        end
    end
end


function constraint_pipe_head_loss(wm::AbstractMILPRModel, n::Int, a::Int, node_fr::Int, node_to::Int, alpha::Float64, L::Float64, r::Float64)
    # If the number of breakpoints is not positive, no constraints are added.
    pipe_breakpoints = get(wm.ext, :pipe_breakpoints, 0)
    if pipe_breakpoints <= 0 return end

    # Add constraints corresponding to positive outer-approximations.
    qp, dhp = var(wm, n, :qp_pipe, a), var(wm, n, :dhp_pipe, a)

    for qp_hat in range(0.0, stop=JuMP.upper_bound(qp), length=pipe_breakpoints)
        lhs = r * _get_head_loss_oa(qp, qp_hat, alpha)
        c_p = JuMP.@constraint(wm.model, lhs <= inv(L) * dhp)
        append!(con(wm, n, :head_loss)[a], [c_p])
    end

    # Add constraints corresponding to positive outer-approximations.
    qn, dhn = var(wm, n, :qn_pipe, a), var(wm, n, :dhn_pipe, a)

    for qn_hat in range(0.0, stop=JuMP.upper_bound(qn), length=pipe_breakpoints)
        lhs = r * _get_head_loss_oa(qn, qn_hat, alpha)
        c_n = JuMP.@constraint(wm.model, lhs <= inv(L) * dhn)
        append!(con(wm, n, :head_loss)[a], [c_n])
    end
end


function constraint_check_valve_head_loss(wm::AbstractMILPRModel, n::Int, a::Int, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    # If the number of breakpoints is not positive, no constraints are added.
    pipe_breakpoints = get(wm.ext, :pipe_breakpoints, 0)
    if pipe_breakpoints <= 0 return end

    # Get common variables for outer approximation constraints.
    qp, z = var(wm, n, :qp_check_valve, a), var(wm, n, :z_check_valve, a)
    dhp = var(wm, n, :dhp_check_valve, a)

    # Add outer approximation constraints.
    for qp_hat in range(0.0, stop=JuMP.upper_bound(qp), length=pipe_breakpoints)
        lhs = _get_head_loss_cv_oa(qp, z, qp_hat, ref(wm, n, :alpha))
        c_p = JuMP.@constraint(wm.model, r * lhs <= inv(L) * dhp)
        append!(con(wm, n, :head_loss)[a], [c_p])
    end
end


function constraint_shutoff_valve_head_loss(wm::AbstractMILPRModel, n::Int, a::Int, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    # If the number of breakpoints is not positive, no constraints are added.
    pipe_breakpoints = get(wm.ext, :pipe_breakpoints, 0)
    if pipe_breakpoints <= 0 return end

    # Get common data for outer approximation constraints.
    qp, qn = var(wm, n, :qp_shutoff_valve, a), var(wm, n, :qn_shutoff_valve, a)
    dhp, dhn = var(wm, n, :dhp_shutoff_valve, a), var(wm, n, :dhn_shutoff_valve, a)
    y = var(wm, n, :y_shutoff_valve, a)

    # Add outer approximation constraints for positively-directed flow.
    for qp_hat in range(0.0, stop=JuMP.upper_bound(qp), length=pipe_breakpoints)
        lhs = _get_head_loss_cv_oa(qp, y, qp_hat, ref(wm, n, :alpha))
        c_p = JuMP.@constraint(wm.model, L * r * lhs <= dhp)
        append!(con(wm, n, :head_loss)[a], [c_p])
    end

    # Add outer approximation constraints for negatively-directed flow.
    for qn_hat in range(0.0, stop=JuMP.upper_bound(qn), length=pipe_breakpoints)
        lhs = _get_head_loss_cv_oa(qn, (1.0 - y), qn_hat, ref(wm, n, :alpha))
        c_n = JuMP.@constraint(wm.model, L * r * lhs <= dhn)
        append!(con(wm, n, :head_loss)[a], [c_n])
    end
end


function objective_owf(wm::AbstractMILPRModel)
    # If the number of breakpoints is not positive, no objective is added.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 0)
    if pump_breakpoints <= 0 return end

    # Initialize the objective function.
    objective = JuMP.AffExpr(0.0)

    for (n, nw_ref) in nws(wm)
        # Get common variables.
        lambda = var(wm, n, :lambda_pump)

        # Get common constant parameters.
        rho = 1000.0 # Water density (kilogram per cubic meter).
        gravity = 9.80665 # Gravitational acceleration (meter per second squared).
        constant = rho * gravity * ref(wm, n, :time_step)

        for (a, pump) in nw_ref[:pump]
            if haskey(pump, "energy_price")
                # Get price and pump curve data.
                price = pump["energy_price"]
                head_curve = ref(wm, n, :pump, a)["head_curve"]
                curve_fun = _get_function_from_head_curve(head_curve)

                # Get flow-related variables and data.
                qp, z = var(wm, n, :qp_pump, a), var(wm, n, :z_pump, a)
                qp_ub = JuMP.upper_bound(qp)

                # Generate a set of uniform flow and cubic function breakpoints.
                breakpoints = range(0.0, stop=qp_ub, length=pump_breakpoints)
                f = _calc_cubic_flow_values(collect(breakpoints), curve_fun)

                # Get pump efficiency data.
                if haskey(pump, "efficiency_curve")
                    eff_curve = pump["efficiency_curve"]
                    eff = _calc_efficiencies(collect(breakpoints), eff_curve)
                else
                    eff = pump["efficiency"]
                end

                # Add the cost corresponding to the current pump's operation.
                inner_expr = (constant*price) .* inv.(eff) .* f
                cost = sum(inner_expr[k]*lambda[a, k] for k in 1:pump_breakpoints)
                JuMP.add_to_expression!(objective, cost)
            else
                Memento.error(_LOGGER, "No cost given for pump \"$(pump["name"])\"")
            end
        end
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end
