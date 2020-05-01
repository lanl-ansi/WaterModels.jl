# Define MILP (approximation-based mixed-integer linear programming)
# implementations of common water distribution model specifications.

function variable_flow_piecewise_weights(wm::AbstractMILPModel; nw::Int=wm.cnw, report::Bool=false)
    # Create weights involved in convex combination constraints for pipes.
    lambda_pipe = var(wm, nw)[:lambda_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe_fixed), k in 1:wm.ext[:pipe_breakpoints]],
        base_name="$(nw)_lambda", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :pipe_fixed, a), "lambda_start", k))

    #report && sol_component_value(wm, nw, :link, :lambda, ids(wm, nw, :pipe_fixed), lambda_pipe)

    # Create weights involved in convex combination constraints for pumps.
    lambda_pump = var(wm, nw)[:lambda_pump] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pump), k in 1:wm.ext[:pump_breakpoints]],
        base_name="$(nw)_lambda", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :pump, a), "lambda_start", k))

    #report && sol_component_value(wm, nw, :link, :lambda, ids(wm, nw, :pump), lambda_pump)
end

function variable_flow_piecewise_weights_des(wm::AbstractMILPModel; nw::Int=wm.cnw, report::Bool=false)
    n_r = Dict(a=>length(ref(wm, nw, :resistance, a)) for a in ids(wm, nw, :link_des))
    K = 1:wm.ext[:pipe_breakpoints]

    # Create weights involved in convex combination constraints.
    lambda_pipe = var(wm, nw)[:lambda_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe_des), r in 1:n_r[a], k in K],
        base_name="$(nw)_lambda", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :pipe_des, a), "lambda_start", r))

    #report && sol_component_value(wm, nw, :link, :lambda, ids(wm, nw, :pipe_des), lambda_pipe)
end

function variable_flow_piecewise_adjacency(wm::AbstractMILPModel; nw::Int=wm.cnw, report::Bool=false)
    # Create binary variables involved in convex combination constraints for pipes.
    x_pw_pipe = var(wm, nw)[:x_pw_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:wm.ext[:pipe_breakpoints]-1],
        base_name="$(nw)_x_pw", binary=true,
        start=comp_start_value(ref(wm, nw, :pipe, a), "x_pw_start"))

    #report && sol_component_value(wm, nw, :link, :x_pw, ids(wm, nw, :pipe), x_pw_pipe)

    # Create binary variables involved in convex combination constraints for pumps.
    x_pw_pump = var(wm, nw)[:x_pw_pump] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pump), k in 1:wm.ext[:pump_breakpoints]-1],
        base_name="$(nw)_x_pw", binary=true,
        start=comp_start_value(ref(wm, nw, :pump, a), "x_pw_start"))

    #report && sol_component_value(wm, nw, :link, :x_pw, ids(wm, nw, :pump), x_pw_pump)
end

"Creates flow variables for `MILP` formulations (`q`, `lambda`, `x_pw`)."
function variable_flow(wm::AbstractMILPModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    # Create common flow variables.
    variable_flow_common(wm, nw=nw, bounded=bounded, report=report)

    # Create variables required for convex combination piecewise approximation.
    variable_flow_piecewise_weights(wm, nw=nw)
    variable_flow_piecewise_adjacency(wm, nw=nw)
end

"Creates network design flow variables for `MILP` formulations (`q_des`, `lambda`, `x_pw`)."
function variable_flow_des(wm::AbstractMILPModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    # Create common flow variables.
    variable_flow_des_common(wm, nw=nw, bounded=bounded, report=report)
    variable_flow_piecewise_weights_des(wm, nw=nw)
end

"Adds head loss constraints for check valves in `MILP` formulations."
function constraint_check_valve_head_loss(wm::AbstractMILPModel, n::Int, a::Int, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    # If the number of breakpoints is not positive, no constraints are added.
    if wm.ext[:pipe_breakpoints] <= 0 return end

    # Gather common variables.
    q, x_cv = [var(wm, n, :q, a), var(wm, n, :x_cv, a)]
    lambda, x_pw = [var(wm, n, :lambda_pipe), var(wm, n, :x_pw_pipe)]
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == x_cv)
    c_2 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == x_cv)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

    # Generate a set of uniform flow and head loss breakpoints.
    q_lb, q_ub = [JuMP.lower_bound(q), JuMP.upper_bound(q)]
    breakpoints = range(q_lb, stop=q_ub, length=wm.ext[:pipe_breakpoints])
    f = _calc_head_loss_values(collect(breakpoints), ref(wm, :alpha))

    # Add a constraint for the flow piecewise approximation.
    K = 1:wm.ext[:pipe_breakpoints]
    q_lhs = sum(breakpoints[k] * lambda[a, k] for k in K)
    c_5 = JuMP.@constraint(wm.model, q_lhs == q)

    # Add a constraint for the head gain piecewise approximation.
    loss = r .* sum(f[k] .* lambda[a, k] for k in K)
    lhs = inv(L) * (h_i - h_j) - loss
    c_6 = JuMP.@constraint(wm.model, lhs <= 0.0)
    dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    c_7 = JuMP.@constraint(wm.model, lhs >= inv(L) * (1.0 - x_cv) * dh_lb)

    # Append the constraint array with the above-generated constraints.
    append!(con(wm, n, :head_loss, a), [c_1, c_2, c_3, c_4, c_5, c_6, c_7])

    # Add the adjacency constraints.
    for k in 2:wm.ext[:pipe_breakpoints]-1
        adjacency = x_pw[a, k-1] + x_pw[a, k]
        c_8_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :head_loss, a), [c_8_k])
    end
end

"Pump head gain constraint when the pump status is ambiguous."
function constraint_pump_head_gain(wm::AbstractMILPModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # If the number of breakpoints is not positive, no constraints are added.
    if wm.ext[:pump_breakpoints] <= 0 return end

    # Gather common variables.
    x_pump = var(wm, n, :x_pump, a)
    q, g = [var(wm, n, :q, a), var(wm, n, :g, a)]
    lambda, x_pw = [var(wm, n, :lambda_pump), var(wm, n, :x_pw_pump)]

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == x_pump)
    c_2 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == x_pump)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

    # Generate a set of uniform flow and head loss breakpoints.
    q_ub = JuMP.upper_bound(q)
    breakpoints = range(0.0, stop=q_ub, length=wm.ext[:pump_breakpoints])
    f = _calc_pump_gain_values(collect(breakpoints), curve_fun)

    # Add a constraint for the flow piecewise approximation.
    K = 1:wm.ext[:pump_breakpoints]
    q_lhs = sum(breakpoints[k] * lambda[a, k] for k in K)
    c_5 = JuMP.@constraint(wm.model, q_lhs == q)

    # Add a constraint for the head gain piecewise approximation.
    lhs = sum(f[k] * lambda[a, k] for k in K)
    c_6 = JuMP.@constraint(wm.model, lhs == g)

    # Append the constraint array with the above-generated constraints.
    append!(con(wm, n, :head_gain, a), [c_1, c_2, c_3, c_4, c_5, c_6])

    # Add the adjacency constraints.
    for k in 2:wm.ext[:pump_breakpoints]-1
        adjacency = x_pw[a, k-1] + x_pw[a, k]
        c_7_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :head_gain, a), [c_7_k])
    end
end

function constraint_head_loss_pipe_des(wm::AbstractMILPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, resistances)
    # If the number of breakpoints is not positive, no constraints are added.
    if wm.ext[:pipe_breakpoints] <= 0 return end

    # Get common variables and data.
    n_b, K = [wm.ext[:pipe_breakpoints], 1:wm.ext[:pipe_breakpoints]]
    lambda, x_pw = [var(wm, n, :lambda_pipe), var(wm, n, :x_pw_pipe)]
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]

    # Initialize flow expression in the head loss relationship.
    q_loss_expr = JuMP.AffExpr(0.0)

    for (r_id, r) in enumerate(resistances)
        q, x_res = [var(wm, n, :q_des, a)[r_id], var(wm, n, :x_res, a)[r_id]]
        q_lb, q_ub = [JuMP.lower_bound(q), JuMP.upper_bound(q)]

        # Add the first required SOS constraints.
        lambda_sum = sum(lambda[a, r_id, k] for k in K)
        c_1 = JuMP.@constraint(wm.model, lambda_sum == x_res)
        c_2 = JuMP.@constraint(wm.model, lambda[a, r_id, 1] <= x_pw[a, 1])
        c_3 = JuMP.@constraint(wm.model, lambda[a, r_id, n_b] <= x_pw[a, n_b-1])

        # Generate a set of uniform flow and head loss breakpoints.
        breakpoints = range(q_lb, stop=q_ub, length=n_b)
        f = _calc_head_loss_values(collect(breakpoints), alpha)

        # Add a constraint for the head loss piecewise approximation.
        expr_r = r .* sum(f[k] .* lambda[a, r_id, k] for k in K)
        JuMP.add_to_expression!(q_loss_expr, expr_r)

        # Add a constraint for the flow piecewise approximation.
        q_lhs = sum(breakpoints[k] * lambda[a, r_id, k] for k in K)
        c_4 = JuMP.@constraint(wm.model, q_lhs == q)

        # Append the constraint array with the above-generated constraints.
        append!(con(wm, n, :head_loss, a), [c_1, c_2, c_3, c_4])

        # Add the adjacency constraints.
        for k in 2:n_b-1
            adjacency = x_pw[a, k-1] + x_pw[a, k]
            c_5_k = JuMP.@constraint(wm.model, lambda[a, r_id, k] <= adjacency)
            append!(con(wm, n, :head_loss, a), [c_5_k])
        end
    end

    # Add the final SOS and approximation constraints.
    c_6 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == 1.0)
    c_7 = JuMP.@constraint(wm.model, q_loss_expr == inv(L) * (h_i - h_j))
    append!(con(wm, n, :head_loss, a), [c_6, c_7])
end

function constraint_head_loss_pipe(wm::AbstractMILPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    # If the number of breakpoints is not positive, no constraints are added.
    if wm.ext[:pipe_breakpoints] <= 0 return end

    # Get required variables.
    q = var(wm, n, :q, a)
    lambda, x_pw = [var(wm, n, :lambda_pipe), var(wm, n, :x_pw_pipe)]
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == 1.0)
    c_2 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == 1.0)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

    # Generate a set of uniform flow and head loss breakpoints.using CPLEX
    q_lb, q_ub = [JuMP.lower_bound(q), JuMP.upper_bound(q)]
    breakpoints = range(q_lb, stop=q_ub, length=wm.ext[:pipe_breakpoints])

    # Add a constraint for the flow piecewise approximation.
    K = 1:wm.ext[:pipe_breakpoints]
    q_lhs = sum(breakpoints[k] * lambda[a, k] for k in K)
    c_5 = JuMP.@constraint(wm.model, q_lhs == q)

    # Add a constraint for the head loss piecewise approximation.
    f = _calc_head_loss_values(collect(breakpoints), alpha)
    lhs = r .* sum(f[k] .* lambda[a, k] for k in K)
    c_6 = JuMP.@constraint(wm.model, lhs == inv(L) * (h_i - h_j))

    # Append the constraint array with the above-generated constraints.
    append!(con(wm, n, :head_loss, a), [c_1, c_2, c_3, c_4, c_5, c_6])

    # Add the adjacency constraints.
    for k in 2:wm.ext[:pipe_breakpoints]-1
        adjacency = x_pw[a, k-1] + x_pw[a, k]
        c_7_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :head_loss, a), [c_7_k])
    end
end

function objective_owf(wm::AbstractMILPModel) 
    # If the number of breakpoints is not positive, no objective is added.
    if wm.ext[:pump_breakpoints] <= 0 return end

    # Initialize the objective function.
    objective = JuMP.AffExpr(0.0)
    K = 1:wm.ext[:pump_breakpoints]
    time_step = wm.ref[:option]["time"]["hydraulic_timestep"]

    for (n, nw_ref) in nws(wm)
        # Get common variables.
        lambda = var(wm, n, :lambda_pump)

        # Get common constant parameters.
        rho = 1000.0 # Water density (kilogram per cubic meter).
        gravity = 9.80665 # Gravitational acceleration (meter per second squared).
        constant = rho * gravity * time_step

        for (a, pump) in nw_ref[:pump]
            if haskey(pump, "energy_price")
                # Get price and pump curve data.
                price = pump["energy_price"]
                pump_curve = ref(wm, n, :pump, a)["pump_curve"]
                curve_fun = _get_function_from_pump_curve(pump_curve)

                # Get flow-related variables and data.
                q, x_pump = [var(wm, n)[:q][a], var(wm, n)[:x_pump][a]]
                q_ub = JuMP.upper_bound(q)

                # Generate a set of uniform flow and cubic function breakpoints.
                breakpoints = range(0.0, stop=q_ub, length=wm.ext[:pump_breakpoints])
                f = _calc_cubic_flow_values(collect(breakpoints), curve_fun)

                # Get pump efficiency data.
                if haskey(pump, "efficiency_curve")
                    eff_curve = pump["efficiency_curve"]
                    eff = _calc_efficiencies(collect(breakpoints), eff_curve)
                else
                    eff = wm.ref[:option]["energy"]["global_efficiency"]
                end

                # Add the cost corresponding to the current pump's operation.
                inner_expr = inv.(eff) .* f
                cost = constant*price*sum(inner_expr[k]*lambda[a, k] for k in K)
                JuMP.add_to_expression!(objective, cost)
            else
                Memento.error(_LOGGER, "No cost given for pump \"$(pump["name"])\"")
            end
        end
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end
