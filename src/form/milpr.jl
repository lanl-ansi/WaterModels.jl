# Define MILP-R (relaxation-based mixed-integer linear programming)
# implementations of common water distribution model specifications.

function _get_owf_oa(q::JuMP.VariableRef, x_pump::JuMP.VariableRef, q_hat::Float64, coeffs::Array{Float64})
    f = coeffs[1]*q_hat^3 + coeffs[2]*q_hat^2 + coeffs[3]*q_hat
    df = 3.0*coeffs[1]*q_hat^2 + 2.0*coeffs[2]*q_hat + coeffs[3]
    return f*x_pump + df*(q - q_hat*x_pump)
end

function _get_head_loss_oa(q::JuMP.VariableRef, q_hat::Float64, alpha::Float64)
    return q_hat^alpha + alpha * q_hat^(alpha - 1.0) * (q - q_hat)
end

function _get_head_loss_cv_oa(q::JuMP.VariableRef, x_cv::JuMP.VariableRef, q_hat::Float64, alpha::Float64)
    return q_hat^alpha*x_cv + alpha * q_hat^(alpha - 1.0) * (q - q_hat*x_cv)
end

function _get_pump_gain_oa(q::JuMP.VariableRef, x_pump::JuMP.VariableRef, q_hat::Float64, curve_fun::Array{Float64})
    f = curve_fun[1]*q_hat^2 + curve_fun[2]*q_hat + curve_fun[3]
    df = 2.0 * curve_fun[1] * q_hat + curve_fun[2]
    return f * x_pump + df * (q - q_hat * x_pump)
end

function variable_flow(wm::AbstractMILPRModel; nw::Int=wm.cnw, bounded::Bool=true)
    # Create directed flow variables (i.e., qp and qn).
    variable_flow_common(wm, nw=nw, bounded=bounded)

    # Variables for direction.
    variable_flow_direction(wm, nw=nw)
end

function variable_pump_operation(wm::AbstractMILPRModel; nw::Int=wm.cnw, report::Bool=true)
    # Create common pump variables.
    variable_pump_common(wm, nw=nw, report=report)

    # Create weights involved in convex combination constraints.
    lambda = var(wm, nw)[:lambda] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pump), k in 1:wm.ext[:num_breakpoints]],
        base_name="$(nw)_lambda", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :pump, a), "lambda_start"))

    # Create binary variables involved in convex combination constraints.
    x_pw = var(wm, nw)[:x_pw] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pump), k in 1:wm.ext[:num_breakpoints]-1],
        base_name="$(nw)_x_pw", binary=true,
        start=comp_start_value(ref(wm, nw, :pump, a), "x_pw_start"))
end

"Pump head gain constraint when the pump status is ambiguous."
function constraint_head_gain_pump(wm::AbstractMILPRModel, n::Int, a::Int, node_fr::Int, node_to::Int, pc::Array{Float64})
    # If the number of breakpoints is not positive, no constraints are added.
    if wm.ext[:num_breakpoints] <= 0 return end

    # Gather flow, head gain, and convex combination variables.
    g = var(wm, n, :g, a)
    qp, x_pump = [var(wm, n, :qp, a), var(wm, n, :x_pump, a)]
    lambda, x_pw = [var(wm, n, :lambda), var(wm, n, :x_pw)]

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == x_pump)
    c_2 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == x_pump)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

    # Add a constraint for the flow piecewise approximation.
    qp_ub, K = [JuMP.upper_bound(qp), 1:wm.ext[:num_breakpoints]]
    breakpoints = range(0.0, stop=qp_ub, length=wm.ext[:num_breakpoints])
    qp_lhs = sum(breakpoints[k] * lambda[a, k] for k in K)
    c_5 = JuMP.@constraint(wm.model, qp_lhs == qp)

    # Add a constraint for the head gain piecewise approximation.
    f = _calc_pump_gain_values(collect(breakpoints), pc)
    lhs = sum(f[k] * lambda[a, k] for k in K)
    c_6 = JuMP.@constraint(wm.model, lhs == g)

    # Append the constraint array with the above-generated constraints.
    append!(con(wm, n, :head_gain, a), [c_1, c_2, c_3, c_4, c_5, c_6])

    # Add the adjacency constraints.
    for k in 2:wm.ext[:num_breakpoints]-1
        adjacency = x_pw[a, k-1] + x_pw[a, k]
        c_7_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :head_gain, a), [c_7_k])
    end
end

function constraint_head_loss_pipe_des(wm::AbstractMILPRModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, resistances)
    # If the number of breakpoints is not positive, no constraints are added.
    if wm.ext[:num_breakpoints] <= 0 return end

    for (r_id, r) in enumerate(resistances)
        # Add constraints corresponding to positive outer-approximations.
        qp, dhp = [var(wm, n, :qp_des, a)[r_id], var(wm, n, :dhp, a)]
        qp_ub = JuMP.upper_bound(qp)

        for qp_hat in range(0.0, stop=qp_ub, length=wm.ext[:num_breakpoints])
            lhs = r * _get_head_loss_oa(qp, qp_hat, alpha)
            c_p = JuMP.@constraint(wm.model, lhs <= inv(L) * dhp)
            append!(con(wm, n, :head_loss)[a], [c_p])
        end

        # Add constraints corresponding to negative outer-approximations.
        qn, dhn = [var(wm, n, :qn_des, a)[r_id], var(wm, n, :dhn, a)]
        qn_ub = JuMP.upper_bound(qn)

        for qn_hat in range(0.0, stop=qn_ub, length=wm.ext[:num_breakpoints])
            lhs = r * _get_head_loss_oa(qn, qn_hat, alpha)
            c_n = JuMP.@constraint(wm.model, lhs <= inv(L) * dhn)
            append!(con(wm, n, :head_loss)[a], [c_n])
        end
    end
end

function constraint_head_loss_pipe(wm::AbstractMILPRModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    # If the number of breakpoints is not positive, no constraints are added.
    if wm.ext[:num_breakpoints] <= 0 return end

    # Add constraints corresponding to positive outer-approximations.
    qp, dhp = [var(wm, n, :qp, a), var(wm, n, :dhp, a)]
    qp_ub = JuMP.upper_bound(qp)

    for qp_hat in range(0.0, stop=qp_ub, length=wm.ext[:num_breakpoints])
        lhs = r * _get_head_loss_oa(qp, qp_hat, alpha)
        c_p = JuMP.@constraint(wm.model, lhs <= inv(L) * dhp)
        append!(con(wm, n, :head_loss)[a], [c_p])
    end

    # Add constraints corresponding to positive outer-approximations.
    qn, dhn = [var(wm, n, :qn, a), var(wm, n, :dhn, a)]
    qn_ub = JuMP.upper_bound(qn)

    for qn_hat in range(0.0, stop=qn_ub, length=wm.ext[:num_breakpoints])
        lhs = r * _get_head_loss_oa(qn, qn_hat, alpha)
        c_n = JuMP.@constraint(wm.model, lhs <= inv(L) * dhn)
        append!(con(wm, n, :head_loss)[a], [c_n])
    end
end

function constraint_head_loss_check_valve(wm::AbstractMILPRModel, n::Int, a::Int, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    # If the number of breakpoints is not positive, no constraints are added.
    if wm.ext[:num_breakpoints] <= 0 return end

    # Get common data for outer-approximation constraints.
    qp, x_cv = [var(wm, n, :qp, a), var(wm, n, :x_cv, a)]
    dhp = var(wm, n, :dhp, a)
    qp_ub = JuMP.upper_bound(qp)

    # Add outer-approximation constraints.
    for qp_hat in range(0.0, stop=qp_ub, length=wm.ext[:num_breakpoints])
        lhs = _get_head_loss_cv_oa(qp, x_cv, qp_hat, ref(wm, n, :alpha))
        c_p = JuMP.@constraint(wm.model, r * lhs <= inv(L) * dhp)
        append!(con(wm, n, :head_loss)[a], [c_p])
    end
end

function objective_owf(wm::AbstractMILPRModel)
    # If the number of breakpoints is not positive, no objective is added.
    if wm.ext[:num_breakpoints] <= 0 return end

    # Initialize the objective function.
    objective = JuMP.AffExpr(0.0)
    K = 1:wm.ext[:num_breakpoints]

    for (n, nw_ref) in nws(wm)
        # Get common variables.
        lambda = var(wm, n, :lambda)

        # Get common constant parameters.
        rho = 1000.0 # Water density (kilogram per cubic meter).
        gravity = 9.80665 # Gravitational acceleration (meter per second squared).
        time_step = nw_ref[:option]["time"]["hydraulic_timestep"]
        constant = rho * gravity * time_step

        for (a, pump) in nw_ref[:pump]
            if haskey(pump, "energy_price")
                # Get price and pump curve data.
                price = pump["energy_price"]
                pump_curve = ref(wm, n, :pump, a)["pump_curve"]
                curve_fun = _get_function_from_pump_curve(pump_curve)

                # Get flow-related variables and data.
                qp, x_pump = [var(wm, n)[:qp][a], var(wm, n)[:x_pump][a]]
                qp_ub = JuMP.upper_bound(qp)

                # Generate a set of uniform flow and cubic function breakpoints.
                breakpoints = range(0.0, stop=qp_ub, length=wm.ext[:num_breakpoints])
                f = _calc_cubic_flow_values(collect(breakpoints), curve_fun)

                # Get pump efficiency data.
                if haskey(pump, "efficiency_curve")
                    eff_curve = pump["efficiency_curve"]
                    eff = _calc_efficiencies(collect(breakpoints), eff_curve)
                else
                    eff = ref(wm, n, :option)["energy"]["global_efficiency"]
                end

                # Add the cost corresponding to the current pump's operation.
                inner_expr = (constant*price) .* inv.(eff) .* f
                cost = sum(inner_expr[k]*lambda[a, k] for k in K)
                JuMP.add_to_expression!(objective, cost)
            else
                Memento.error(_LOGGER, "No cost given for pump \"$(pump["name"])\"")
            end
        end
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end
