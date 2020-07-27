# Define common MICP (mixed-integer convex program) implementations of water
# distribution constraints, which use directed flow variables.

function variable_pump_operation(wm::AbstractMICPModel; nw::Int=wm.cnw, report::Bool=true)
    # If the number of breakpoints is zero, the variables below are not added.
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

function constraint_pipe_head_loss_des(wm::AbstractMICPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, pipe_resistances)
    # Collect head difference variables.
    dhp, dhn = var(wm, n, :dhp_pipe, a), var(wm, n, :dhn_pipe, a)

    for (r_id, r) in enumerate(pipe_resistances)
        # Collect directed flow variables.
        qp, qn = [var(wm, n, :qp_des, a)[r_id], var(wm, n, :qn_des, a)[r_id]]

        # Build the relaxed head loss constraints.
        c_p = JuMP.@NLconstraint(wm.model, r * head_loss(qp) <= inv(L) * dhp)
        c_n = JuMP.@NLconstraint(wm.model, r * head_loss(qn) <= inv(L) * dhn)

        # Append the constraint array.
        append!(con(wm, n, :head_loss, a), [c_p, c_n])
    end
end

"Pump head gain constraint when the pump status is ambiguous."
function constraint_pump_head_gain(wm::AbstractMICPModel, n::Int, a::Int, node_fr::Int, node_to::Int, pc::Array{Float64})
    # Gather flow and head gain variables.
    g = var(wm, n, :g, a)
    qp, z = var(wm, n, :qp_pump, a), var(wm, n, :z_pump, a)

    # Gather head-related variables and data.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)
    dhp, dhn = var(wm, n, :dhp_pump, a), var(wm, n, :dhn_pump, a)
    dhp_ub, dhn_ub = JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)

    # Define the (relaxed) head gain caused by the pump.
    g_expr = pc[1]*qp^2 + pc[2]*qp + pc[3]*z
    c = JuMP.@constraint(wm.model, g_expr >= g) # Concavified.

    # Append the constraint array.
    append!(con(wm, n, :head_gain, a), [c])

    # If the number of breakpoints is not positive, no constraints are added.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 0)
    if pump_breakpoints <= 0 return end

    # Gather flow, head gain, and convex combination variables.
    lambda, x_pw = [var(wm, n, :lambda_pump), var(wm, n, :x_pw_pump)]

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == z)
    c_2 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == z)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

    # Add a constraint for the flow piecewise approximation.
    qp_ub = JuMP.upper_bound(qp)
    breakpoints = range(0.0, stop=qp_ub, length=pump_breakpoints)
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
end

function constraint_check_valve_head_loss(wm::AbstractMICPModel, n::Int, a::Int, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    # Gather flow- and check valve-related variables.
    qp, qn = var(wm, n, :qp_check_valve, a), var(wm, n, :qn_check_valve, a)
    z = var(wm, n, :z_check_valve, a)

    # Gather head variables and upper bound data.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)
    dhp, dhn = var(wm, n, :dhp_check_valve, a), var(wm, n, :dhn_check_valve, a)
    dhp_ub, dhn_ub = JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)

    # Add constraints for flow in the positive and negative directions.
    lhs = JuMP.@NLexpression(wm.model, r*head_loss(qp) - inv(L)*dhp)
    c_p = JuMP.@NLconstraint(wm.model, lhs <= inv(L) * dhp_ub * (1.0 - z))
    c_n = JuMP.@NLconstraint(wm.model, dhn <= inv(L) * dhn_ub * (1.0 - z))

    # Append the constraint array.
    append!(con(wm, n, :head_loss)[a], [c_p, c_n])
end

function constraint_shutoff_valve_head_loss(wm::AbstractMICPModel, n::Int, a::Int, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    # Gather flow- and shutoff valve-related variables.
    qp, qn = var(wm, n, :qp_shutoff_valve, a), var(wm, n, :qn_shutoff_valve, a)
    y, z = var(wm, n, :y_shutoff_valve, a), var(wm, n, :z_shutoff_valve, a)

    # Gather head variables and upper bound data.
    dhp, dhn = var(wm, n, :dhp_shutoff_valve, a), var(wm, n, :dhn_shutoff_valve, a)
    dhp_ub, dhn_ub = JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)

    # Add constraints for flow in the positive and negative directions.
    lhs_p = JuMP.@NLexpression(wm.model, L * r * head_loss(qp) - dhp)
    c_1 = JuMP.@NLconstraint(wm.model, lhs_p <= dhp_ub * (1.0 - y))
    c_2 = JuMP.@NLconstraint(wm.model, lhs_p <= dhp_ub * (1.0 - z))

    lhs_n = JuMP.@NLexpression(wm.model, L * r * head_loss(qn) - dhn)
    c_3 = JuMP.@NLconstraint(wm.model, lhs_n <= dhn_ub * y)
    c_4 = JuMP.@NLconstraint(wm.model, lhs_n <= dhn_ub * (1.0 - z))

    # Append the constraint array.
    append!(con(wm, n, :head_loss)[a], [c_1, c_2, c_3, c_4])
end

function constraint_pipe_head_loss(wm::AbstractMICPModel, n::Int, a::Int, node_fr::Int, node_to::Int, alpha::Float64, L::Float64, r::Float64)
    # Gather common variables.
    qp, qn = var(wm, n, :qp_pipe, a), var(wm, n, :qn_pipe, a)
    dhp, dhn = var(wm, n, :dhp_pipe, a), var(wm, n, :dhn_pipe, a)

    # Add constraints for head loss in the positive and negative directions.
    c_p = JuMP.@NLconstraint(wm.model, r * head_loss(qp) <= inv(L) * dhp)
    c_n = JuMP.@NLconstraint(wm.model, r * head_loss(qn) <= inv(L) * dhn)

    # Append the constraint array.
    append!(con(wm, n, :head_loss)[a], [c_p, c_n])
end

function objective_wf(wm::AbstractMICPModel)
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end

function objective_owf(wm::AbstractMICPModel)
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
