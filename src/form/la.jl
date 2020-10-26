# Define LA (linear approximation-based) implementations of water models.

function variable_flow_piecewise_weights(wm::LAWaterModel; nw::Int=wm.cnw, report::Bool=false)
    # Get the number of breakpoints for the pipe.
    pipe_breakpoints = get(wm.ext, :pipe_breakpoints, 2)

    # Create weights involved in convex combination constraints for pipes.
    lambda_pipe = var(wm, nw)[:lambda_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:pipe_breakpoints],
        base_name="$(nw)_lambda", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :pipe, a), "lambda_start", k))

    # Create weights involved in convex combination constraints.
    n_r = Dict(a=>length(ref(wm, nw, :resistance, a)) for a in ids(wm, nw, :des_pipe))
    lambda_des_pipe = var(wm, nw)[:lambda_des_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :des_pipe), r in 1:n_r[a], k in 1:pipe_breakpoints],
        base_name="$(nw)_lambda", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :des_pipe, a), "lambda_start", r))

    # Get the number of breakpoints for the pump.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 2)

    # Create weights involved in convex combination constraints for pumps.
    lambda_pump = var(wm, nw)[:lambda_pump] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pump), k in 1:pump_breakpoints],
        base_name="$(nw)_lambda", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :pump, a), "lambda_start", k))
end


function variable_flow_piecewise_adjacency(wm::LAWaterModel; nw::Int=wm.cnw, report::Bool=false)
    # Get the number of breakpoints for the pipe.
    pipe_breakpoints = get(wm.ext, :pipe_breakpoints, 2)

    # Create binary variables for pipe convex combination constraints.
    x_pw_pipe = var(wm, nw)[:x_pw_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:pipe_breakpoints-1], base_name="$(nw)_x_pw", binary=true,
        start=comp_start_value(ref(wm, nw, :pipe, a), "x_pw_start"))

    # Create binary variables for design pipe convex combination constraints.
    x_pw_des_pipe = var(wm, nw)[:x_pw_des_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :des_pipe), k in 1:pipe_breakpoints-1], base_name="$(nw)_x_pw", binary=true,
        start=comp_start_value(ref(wm, nw, :des_pipe, a), "x_pw_start"))

    # Get the number of breakpoints for the pump.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 2)

    # Create binary variables involved in convex combination constraints for pumps.
    x_pw_pump = var(wm, nw)[:x_pw_pump] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pump), k in 1:pump_breakpoints-1], base_name="$(nw)_x_pw",
        binary=true, start=comp_start_value(ref(wm, nw, :pump, a), "x_pw_start"))
end


"Creates flow variables for `LA` formulations (`q`, `lambda`, `x_pw`)."
function variable_flow(wm::LAWaterModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    for name in ["pipe", "pump", "regulator", "short_pipe", "valve"]
        # Create directed flow (`qp` and `qn`) variables for each component.
        variable_component_flow(wm, name; nw=nw, bounded=bounded, report=report)
    end

    # Create variables required for convex combination piecewise approximation.
    variable_flow_piecewise_weights(wm, nw=nw)
    variable_flow_piecewise_adjacency(wm, nw=nw)

    # Create common flow variables.
    variable_flow_des_common(wm, nw=nw, bounded=bounded, report=report)
end


"Pump head gain constraint when the pump status is ambiguous."
function constraint_on_off_pump_head_gain(wm::LAWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64}, q_min_forward::Float64)
    # Get the number of breakpoints for the pump.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 2)

    # Gather common variables.
    q, g, z = var(wm, n, :q_pump, a), var(wm, n, :g_pump, a), var(wm, n, :z_pump, a)
    lambda, x_pw = var(wm, n, :lambda_pump), var(wm, n, :x_pw_pump)

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == z)
    c_2 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == z)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

    # Generate a set of uniform flow breakpoints.
    breakpoints = range(q_min_forward, stop=JuMP.upper_bound(q), length=pump_breakpoints)

    # Add a constraint for the flow piecewise approximation.
    q_lhs = sum(breakpoints[k] * lambda[a, k] for k in 1:pump_breakpoints)
    c_5 = JuMP.@constraint(wm.model, q_lhs == q)

    # Add a constraint for the head gain piecewise approximation.
    f = _calc_pump_gain_values(collect(breakpoints), curve_fun)
    g_lhs = sum(f[k] * lambda[a, k] for k in 1:pump_breakpoints)
    c_6 = JuMP.@constraint(wm.model, g_lhs == g)

    # Append the constraint array with the above-generated constraints.
    append!(con(wm, n, :on_off_pump_head_gain, a), [c_1, c_2, c_3, c_4, c_5, c_6])

    for k in 2:pump_breakpoints-1
        # Add adjacency constraints for each interval.
        adjacency = x_pw[a, k-1] + x_pw[a, k]
        c_7_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :on_off_pump_head_gain, a), [c_7_k])
    end
end


function constraint_pipe_head_loss(wm::LAWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64, L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # If the number of breakpoints is not positive, no constraints are added.
    pipe_breakpoints = get(wm.ext, :pipe_breakpoints, 2)

    # Get required variables.
    q, h_i, h_j = var(wm, n, :q_pipe, a), var(wm, n, :h, node_fr), var(wm, n, :h, node_to)
    lambda, x_pw = var(wm, n, :lambda_pipe), var(wm, n, :x_pw_pipe)

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == 1.0)
    c_2 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == 1.0)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

    # Generate a set of uniform flow breakpoints.
    q_lb, q_ub = JuMP.lower_bound(q), JuMP.upper_bound(q)
    breakpoints = range(q_lb, stop=q_ub, length=pipe_breakpoints)

    # Add a constraint for the flow piecewise approximation.
    q_lhs = sum(breakpoints[k] * lambda[a, k] for k in 1:pipe_breakpoints)
    c_5 = JuMP.@constraint(wm.model, q_lhs == q)

    # Add a constraint for the head loss piecewise approximation.
    f = _calc_head_loss_values(collect(breakpoints), exponent)
    lhs = r * sum(f[k] * lambda[a, k] for k in 1:pipe_breakpoints)
    c_6 = JuMP.@constraint(wm.model, lhs == inv(L) * (h_i - h_j))

    # Append the constraint array with the above-generated constraints.
    append!(con(wm, n, :pipe_head_loss, a), [c_1, c_2, c_3, c_4, c_5, c_6])

    # Add the adjacency constraints.
    for k in 2:pipe_breakpoints-1
        adjacency = x_pw[a, k-1] + x_pw[a, k]
        c_7_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :pipe_head_loss, a), [c_7_k])
    end
end


"Adds head loss constraint for a design pipe in the `LA` formulation."
function constraint_on_off_pipe_head_loss_des(wm::LAWaterModel, n::Int, a::Int, exponent::Float64, node_fr::Int, node_to::Int, L::Float64, resistances)
    # Get the number of breakpoints for the pipe.
    n_b = pipe_breakpoints = get(wm.ext, :pipe_breakpoints, 2)

    # Get common variables and data.
    lambda, x_pw = var(wm, n, :lambda_des_pipe), var(wm, n, :x_pw_des_pipe)
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # Initialize flow expression in the head loss relationship.
    q_loss_expression = JuMP.AffExpr(0.0)

    for (r_id, r) in enumerate(resistances)
        q, z = var(wm, n, :q_des_pipe, a)[r_id], var(wm, n, :z_des_pipe, a)[r_id]
        q_lb, q_ub = JuMP.lower_bound(q), JuMP.upper_bound(q)

        # Add the first required SOS constraints.
        lambda_sum = sum(lambda[a, r_id, k] for k in 1:pipe_breakpoints)
        c_1 = JuMP.@constraint(wm.model, lambda_sum == z)
        c_2 = JuMP.@constraint(wm.model, lambda[a, r_id, 1] <= x_pw[a, 1])
        c_3 = JuMP.@constraint(wm.model, lambda[a, r_id, n_b] <= x_pw[a, n_b-1])

        # Generate a set of uniform flow breakpoints.
        breakpoints = range(q_lb, stop=q_ub, length=n_b)

        # Add a constraint for the head loss piecewise approximation.
        f = _calc_head_loss_values(collect(breakpoints), exponent)
        expr_r = r .* sum(f[k] .* lambda[a, r_id, k] for k in 1:pipe_breakpoints)
        JuMP.add_to_expression!(q_loss_expression, expr_r)

        # Add a constraint for the flow piecewise approximation.
        q_lhs = sum(breakpoints[k] * lambda[a, r_id, k] for k in 1:pipe_breakpoints)
        c_4 = JuMP.@constraint(wm.model, q_lhs == q)

        # Append the constraint array with the above-generated constraints.
        append!(con(wm, n, :on_off_pipe_head_loss_des, a), [c_1, c_2, c_3, c_4])

        # Add the adjacency constraints.
        for k in 2:n_b-1
            adjacency = x_pw[a, k-1] + x_pw[a, k]
            c_5_k = JuMP.@constraint(wm.model, lambda[a, r_id, k] <= adjacency)
            append!(con(wm, n, :on_off_pipe_head_loss_des, a), [c_5_k])
        end
    end

    # Add the final SOS and approximation constraints.
    c_6 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == 1.0)
    c_7 = JuMP.@constraint(wm.model, q_loss_expression == inv(L) * (h_i - h_j))
    append!(con(wm, n, :on_off_pipe_head_loss_des, a), [c_6, c_7])
end


function objective_owf(wm::LAWaterModel) 
    # Get the number of breakpoints for the pump.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 2)

    # Initialize the objective function.
    objective = JuMP.AffExpr(0.0)

    for (n, nw_ref) in nws(wm)
        # Get common variables.
        lambda = var(wm, n, :lambda_pump)

        # Get common constant parameters.
        constant = _DENSITY * _GRAVITY * ref(wm, n, :time_step)

        for (a, pump) in nw_ref[:pump]
            if haskey(pump, "energy_price")
                # Get price and pump curve data.
                price = pump["energy_price"]
                head_curve = ref(wm, n, :pump, a)["head_curve"]
                curve_fun = _get_function_from_head_curve(head_curve)

                # Get flow-related variables and data.
                q, z = var(wm, n, :q_pump, a), var(wm, n, :z_pump, a)
                q_ub = JuMP.upper_bound(q)

                # Generate a set of uniform flow and cubic function breakpoints.
                breakpoints = range(0.0, stop=q_ub, length=pump_breakpoints)
                f = _calc_cubic_flow_values(collect(breakpoints), curve_fun)

                # Get pump efficiency data.
                if haskey(pump, "efficiency_curve")
                    eff_curve = pump["efficiency_curve"]
                    eff = _calc_efficiencies(collect(breakpoints), eff_curve)
                else
                    eff = pump["efficiency"]
                end

                # Add the cost corresponding to the current pump's operation.
                inner_expr = constant * price * inv.(eff) .* f
                cost = sum(inner_expr[k]*lambda[a, k] for k in 1:pump_breakpoints)
                JuMP.add_to_expression!(objective, cost)
            else
                Memento.error(_LOGGER, "No cost given for pump \"$(pump["name"])\"")
            end
        end
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end
