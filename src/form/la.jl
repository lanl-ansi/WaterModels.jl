# Define LA (linear approximation-based) implementations of water models.

function variable_flow_piecewise_weights(wm::LAWaterModel; nw::Int=wm.cnw, report::Bool=false)
    # Get the number of breakpoints for the pipe.
    pipe_breakpoints = get(wm.ext, :pipe_breakpoints, 2)

    # Create weights involved in convex combination constraints for pipes.
    lambda_pipe = var(wm, nw)[:lambda_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pipe), k in 1:pipe_breakpoints],
        base_name="$(nw)_lambda", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :pipe, a), "lambda_start", k))

    # Create weights involved in convex combination constraints for design pipes.
    lambda_des_pipe = var(wm, nw)[:lambda_des_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :des_pipe), k in 1:pipe_breakpoints],
        base_name="$(nw)_lambda", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :des_pipe, a), "lambda_start", k))

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
    for name in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
        # Create flow variables for each node-connecting component.
        _variable_component_flow(wm, name; nw=nw, bounded=bounded, report=report)
    end

    # Create variables required for convex combination piecewise approximation.
    variable_flow_piecewise_weights(wm; nw = nw)
    variable_flow_piecewise_adjacency(wm; nw = nw)
end


function constraint_pipe_head_loss(
    wm::LAWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
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
    breakpoints = range(q_lb, stop = q_ub, length = pipe_breakpoints)

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


function constraint_on_off_pump_head_gain(wm::LAWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, q_min_forward::Float64)
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
    head_curve_function = _calc_head_curve_function(ref(wm, n, :pump, a))

    # Add a constraint that linearly approximates the head gain variable.
    f = head_curve_function.(collect(breakpoints))
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


function constraint_on_off_des_pipe_head_loss(wm::LAWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64, L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # If the number of breakpoints is not positive, no constraints are added.
    pipe_breakpoints = get(wm.ext, :pipe_breakpoints, 2)

    # Get required variables.
    q, z = var(wm, n, :q_des_pipe, a), var(wm, n, :z_des_pipe, a)
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)
    lambda, x_pw = var(wm, n, :lambda_des_pipe), var(wm, n, :x_pw_des_pipe)

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == z)
    c_2 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == z)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

    # Generate a set of uniform flow breakpoints.
    q_lb, q_ub = JuMP.lower_bound(q), JuMP.upper_bound(q)
    breakpoints = range(q_lb, stop = q_ub, length = pipe_breakpoints)

    # Add a constraint for the flow piecewise approximation.
    q_lhs = sum(breakpoints[k] * lambda[a, k] for k in 1:pipe_breakpoints)
    c_5 = JuMP.@constraint(wm.model, q_lhs == q)

    # Add a constraint for the head loss piecewise approximation.
    f = _calc_head_loss_values(collect(breakpoints), exponent)
    lhs = r * sum(f[k] * lambda[a, k] for k in 1:pipe_breakpoints)

    # TODO: Use a McCormick expansion of the below multiplication with z.
    c_6 = JuMP.@constraint(wm.model, lhs == inv(L) * (h_i - h_j) * z)

    # Append the constraint array with the above-generated constraints.
    append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_1, c_2, c_3, c_4, c_5, c_6])

    # Add the adjacency constraints.
    for k in 2:pipe_breakpoints-1
        adjacency = x_pw[a, k-1] + x_pw[a, k]
        c_7_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :on_off_des_pipe_head_loss, a), [c_7_k])
    end
end


function objective_owf_default(wm::LAWaterModel) 
    # Get the number of breakpoints for the pump.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 2)

    # Initialize the objective function.
    objective = JuMP.AffExpr(0.0)

    for (n, nw_ref) in nws(wm)
        # Get convex combination variable set.
        lambda = var(wm, n, :lambda_pump)

        for (a, pump) in nw_ref[:pump]
            # Ensure that the pump has an associated energy price.
            @assert haskey(pump, "energy_price")

            # Add piecewise linear approximation of cost to the objective.
            q_bp, f_bp = _calc_pump_energy_points(wm, n, a, pump_breakpoints)
            energy = sum(f_bp[k] * lambda[a, k] for k in 1:pump_breakpoints)
            JuMP.add_to_expression!(objective, pump["energy_price"] * energy)
        end
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end
