# Define MILP (approximation-based mixed-integer linear programming)
# implementations of water distribution models.

function variable_flow_piecewise_weights(wm::AbstractMILPModel; nw::Int=wm.cnw, report::Bool=false)
    lambda = var(wm, nw)[:lambda] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :link_fixed), k in 1:wm.ext[:num_breakpoints]],
        base_name="$(nw)_lambda", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :link_fixed, a), "lambda_start"))

    # Initialize the solution reporting data structures.
    report && sol_component_value(wm, nw, :link, :lambda, ids(wm, nw, :link_fixed), lambda)
end

function variable_flow_piecewise_weights_ne(wm::AbstractMILPModel, n::Int=wm.cnw)
    # Set the number of breakpoints used in each outer-approximation.
    n_r = Dict(a=>length(ref(wm, nw, :resistance, a)) for a in ids(wm, nw, :link_ne))
    K = 1:wm.ext[:num_breakpoints]

    lambda = var(wm, nw)[:lambda] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :link_ne), r in 1:n_r[a], k in K],
        base_name="lambda[$(n)]", lower_bound=0.0, upper_bound=1.0,
        start=get_start(ref(wm, nw, :link), a, "lambda_start"))

    # Initialize the solution reporting data structures.
    report && sol_component_value(wm, nw, :link, :lambda, ids(wm, nw, :link_fixed), lambda)
end

function variable_flow_piecewise_adjacency(wm::AbstractMILPModel; nw::Int=wm.cnw, report::Bool=false)
    # Set the number of breakpoints used in each outer-approximation.
    x_pw = var(wm, nw)[:x_pw] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :link_fixed), k in 1:wm.ext[:num_breakpoints]-1],
        base_name="$(nw)_x_pw", binary=true,
        start=comp_start_value(ref(wm, nw, :link_fixed, a), "x_pw_start"))

    # Initialize the solution reporting data structures.
    report && sol_component_value(wm, nw, :link, :x_pw, ids(wm, nw, :link_fixed), x_pw)
end

function variable_flow_piecewise_adjacency_ne(wm::AbstractMILPModel, n::Int=wm.cnw)
    # Set the number of breakpoints used in each outer-approximation.
    var(wm, nw)[:x_pw] = JuMP.@variable(wm.model, [a in ids(wm, nw, :link_ne),
        k in 1:wm.ext[:num_breakpoints]-1], base_name="x_pw[$(n)]", binary=true,
        start=get_start(ref(wm, nw, :link_ne), a, "x_pw_start", 0.0))
end

function variable_flow(wm::AbstractMILPModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    # Create common flow variables.
    variable_flow_common(wm, nw=nw, bounded=bounded, report=report)

    # Create variables required for convex combination piecewise approximation.
    variable_flow_piecewise_weights(wm, nw=nw)
    variable_flow_piecewise_adjacency(wm, nw=nw)
end

"Create network expansion flow variables for undirected flow formulations."
function variable_flow_ne(wm::AbstractMILPModel; nw::Int=wm.cnw, bounded::Bool=true)
    # Create undirected flow variables (i.e., q).
    bounded ? variable_flow_bounded_ne(wm, nw=nw) :
        variable_flow_unbounded_ne(wm, nw=nw)

    # Create expressions capturing the relationships among q, and q_ne.
    var(wm, nw)[:q] = JuMP.@expression(wm.model, [a in ids(wm, nw, :link_ne)],
        sum(var(wm, nw, :q_ne, a)))

    # Create resistance binary variables.
    variable_resistance(wm, nw=nw)

    # Create variables required for convex combination piecewise approximation.
    variable_flow_piecewise_weights_ne(wm, nw)
    variable_flow_piecewise_adjacency_ne(wm, nw)
end

function _get_pump_gain_values(breakpoints::Array{Float64}, curve_fun::Array{Float64})
    return [curve_fun[1]*x*x + curve_fun[2]*x + curve_fun[3] for x in breakpoints]
end

"Pump head gain constraint when the pump status is ambiguous."
function constraint_head_gain_pump(wm::AbstractMILPModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Gather common variables.
    x_pump = var(wm, n, :x_pump, a)
    q, g = [var(wm, n, :q, a), var(wm, n, :g, a)]
    lambda, x_pw = [var(wm, n, :lambda), var(wm, n, :x_pw)]
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]

    # If the pump is off, the flow along the pump must be zero.
    q_lb, q_ub = [JuMP.lower_bound(q), JuMP.upper_bound(q)]
    c_1 = JuMP.@constraint(wm.model, q <= q_ub * x_pump)
    c_2 = JuMP.@constraint(wm.model, q >= 6.31465679e-6 * x_pump)

    # If the pump is off, decouple the head difference relationship.
    dh_lb = JuMP.lower_bound(h_j) - JuMP.upper_bound(h_i)
    c_3 = JuMP.@constraint(wm.model, h_j - h_i - g >= dh_lb * (1.0 - x_pump))
    dh_ub = JuMP.upper_bound(h_j) - JuMP.lower_bound(h_i)
    c_4 = JuMP.@constraint(wm.model, h_j - h_i - g <= dh_ub * (1.0 - x_pump))

    # Append the constraint array.
    append!(con(wm, n, :head_gain, a), [c_1, c_2, c_3, c_4])

    # If the number of breakpoints is not positive, no constraints are added.
    if wm.ext[:num_breakpoints] <= 0 return end

    # Add the required SOS constraints.
    c_5 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == x_pump)
    c_6 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == x_pump)
    c_7 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
    c_8 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

    # Generate a set of uniform flow and head loss breakpoints.
    breakpoints = range(q_lb, stop=q_ub, length=wm.ext[:num_breakpoints])
    f = _get_pump_gain_values(collect(breakpoints), curve_fun)

    # Add a constraint for the flow piecewise approximation.
    K = 1:wm.ext[:num_breakpoints]
    q_lhs = sum(breakpoints[k] * lambda[a, k] for k in K)
    c_9 = JuMP.@constraint(wm.model, q_lhs == q)

    # Add a constraint for the head gain piecewise approximation.
    lhs = sum(f[k] * lambda[a, k] for k in 1:wm.ext[:num_breakpoints])
    c_10 = JuMP.@constraint(wm.model, lhs == g)

    # Append the constraint array with the above-generated constraints.
    append!(con(wm, n, :head_gain, a), [c_5, c_6, c_7, c_8, c_9, c_10])

    # Add the adjacency constraints.
    for k in 2:wm.ext[:num_breakpoints]-1
        adjacency = x_pw[a, k-1] + x_pw[a, k]
        c_11_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :head_gain, a), [c_11_k])
    end
end

function constraint_head_loss_pipe_ne(wm::AbstractMILPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, resistances)
    # Set the number of breakpoints used in each outer-approximation.
    n_b = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 2

    if n_b > 0
        h_i = var(wm, n, :h, node_fr)
        h_j = var(wm, n, :h, node_to)

        lambda = var(wm, n, :lambda)
        x_pw = var(wm, n, :x_pw)
        lhs = JuMP.AffExpr(0.0)

        for (r_id, r) in enumerate(resistances)
            x_res = var(wm, n, :x_res, a)[r_id]
            q_ne = var(wm, n, :q_ne, a)[r_id]

            c_1 = JuMP.@constraint(wm.model, sum(lambda[a, r_id, k] for k in 1:n_b) == x_res)
            c_2 = JuMP.@constraint(wm.model, lambda[a, r_id, 1] <= x_pw[a, 1])
            c_3 = JuMP.@constraint(wm.model, lambda[a, r_id, n_b] <= x_pw[a, n_b-1])
            append!(con(wm, n, :head_loss, a), [c_1, c_2, c_3])

            for k in 2:n_b-1
                c_4 = JuMP.@constraint(wm.model, lambda[a, r_id, k] <= x_pw[a, k-1] + x_pw[a, k])
                append!(con(wm, n, :head_loss, a), [c_4])
            end

            q_ne = var(wm, n, :q_ne, a)[r_id]
            q_ne_lb = JuMP.has_lower_bound(q_ne) ? JuMP.lower_bound(q_ne) : -10.0
            q_ne_ub = JuMP.has_upper_bound(q_ne) ? JuMP.upper_bound(q_ne) : 10.0

            breakpoints = range(q_ne_lb, stop=q_ne_ub, length=n_b)
            f = get_head_loss_values(collect(breakpoints), alpha)
            JuMP.add_to_expression!(lhs, r * sum(f[k] * lambda[a, r_id, k] for k in 1:n_b))

            q_ne_lhs = sum(breakpoints[k] * lambda[a, r_id, k] for k in 1:n_b)
            c_5 = JuMP.@constraint(wm.model, q_ne_lhs == q_ne)
            append!(con(wm, n, :head_loss, a), [c_5])
        end

        c_6 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == 1.0)
        c_7 = JuMP.@constraint(wm.model, lhs == inv(L) * (h_i - h_j))
        append!(con(wm, n, :head_loss, a), [c_6, c_7])
    end
end

function _get_head_loss_values(breakpoints::Array{Float64}, alpha::Float64)
    return [sign(x) * abs(x)^alpha for x in breakpoints]
end

function constraint_head_loss_check_valve(wm::AbstractMILPModel, n::Int, a::Int, node_fr::Int, node_to::Int, L::Float64, r::Float64) 
    # Gather common variables.
    q, x_cv = [var(wm, n, :q, a), var(wm, n, :x_cv, a)]
    lambda, x_pw = [var(wm, n, :lambda), var(wm, n, :x_pw)]
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]

    # If the pump is off, the flow along the pump must be zero.
    q_lb, q_ub = [JuMP.lower_bound(q), JuMP.upper_bound(q)]
    c_1 = JuMP.@constraint(wm.model, q <= q_ub * x_cv)
    c_2 = JuMP.@constraint(wm.model, q >= 6.31465679e-6 * x_cv)

    # Append the constraint array.
    append!(con(wm, n, :head_loss, a), [c_1, c_2])

    # If the number of breakpoints is not positive, no constraints are added.
    if wm.ext[:num_breakpoints] <= 0 return end

    # Add the required SOS constraints.
    c_3 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == x_cv)
    c_4 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == x_cv)
    c_5 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
    c_6 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

    # Generate a set of uniform flow and head loss breakpoints.
    breakpoints = range(q_lb, stop=q_ub, length=wm.ext[:num_breakpoints])
    f = _get_head_loss_values(collect(breakpoints), ref(wm, :alpha))

    # Add a constraint for the flow piecewise approximation.
    K = 1:wm.ext[:num_breakpoints]
    q_lhs = sum(breakpoints[k] * lambda[a, k] for k in K)
    c_7 = JuMP.@constraint(wm.model, q_lhs == q)

    # Add a constraint for the head gain piecewise approximation.
    loss = r .* sum(f[k] .* lambda[a, k] for k in 1:wm.ext[:num_breakpoints])
    lhs = inv(L) * (h_i - h_j) - loss
    c_8 = JuMP.@constraint(wm.model, lhs <= 0.0)
    dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    c_9 = JuMP.@constraint(wm.model, lhs >= inv(L) * (1.0 - x_cv) * dh_lb)

    # Append the constraint array with the above-generated constraints.
    append!(con(wm, n, :head_loss, a), [c_3, c_4, c_5, c_6, c_7, c_8, c_9])

    # Add the adjacency constraints.
    for k in 2:wm.ext[:num_breakpoints]-1
        adjacency = x_pw[a, k-1] + x_pw[a, k]
        c_10_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :head_loss, a), [c_10_k])
    end
end

function constraint_head_loss_pipe(wm::AbstractMILPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    # If the number of breakpoints is not positive, no constraints are added.
    if wm.ext[:num_breakpoints] <= 0 return end

    # Get required variables.
    q = var(wm, n, :q, a)
    lambda, x_pw = [var(wm, n, :lambda), var(wm, n, :x_pw)]
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]

    # Add the required SOS constraints.
    c_1 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == 1.0)
    c_2 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == 1.0)
    c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
    c_4 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

    # Generate a set of uniform flow and head loss breakpoints.
    q_lb, q_ub = [JuMP.lower_bound(q), JuMP.upper_bound(q)]
    breakpoints = range(q_lb, stop=q_ub, length=wm.ext[:num_breakpoints])
    f = _get_head_loss_values(collect(breakpoints), alpha)

    # Add a constraint for the head loss piecewise approximation.
    lhs = r .* sum(f[k] .* lambda[a, k] for k in 1:wm.ext[:num_breakpoints])
    c_5 = JuMP.@constraint(wm.model, lhs == inv(L) * (h_i - h_j))

    # Add a constraint for the flow piecewise approximation.
    K = 1:wm.ext[:num_breakpoints]
    q_lhs = sum(breakpoints[k] * lambda[a, k] for k in K)
    c_6 = JuMP.@constraint(wm.model, q_lhs == q)

    # Append the constraint array with the above-generated constraints.
    append!(con(wm, n, :head_loss, a), [c_1, c_2, c_3, c_4, c_5, c_6])

    # Add the adjacency constraints.
    for k in 2:wm.ext[:num_breakpoints]-1
        adjacency = x_pw[a, k-1] + x_pw[a, k]
        c_7_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :head_loss, a), [c_7_k])
    end
end

function objective_wf(wm::AbstractMILPModel)
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end

function _get_cubic_flow_values(breakpoints::Array{Float64}, curve_fun::Array{Float64})
    return [curve_fun[1]*x^3 + curve_fun[2]*x^2 + curve_fun[3]*x for x in breakpoints]
end

function objective_owf(wm::AbstractMILPModel) 
    # If the number of breakpoints is not positive, no objective is added.
    if wm.ext[:num_breakpoints] <= 0 return end

    # Initialize the objective function.
    objective = JuMP.AffExpr(0.0)

    for (n, nw_ref) in nws(wm)
        # Get common variables.
        lambda = var(wm, n, :lambda)

        # Get common constant parameters.
        efficiency = 1.00 # TODO: Change this after discussion. 0.85 follows Fooladivanda.
        rho = 1000.0 # Water density (kilogram per cubic meter).
        gravity = 9.80665 # Gravitational acceleration (meter per second squared).
        time_step = nw_ref[:option]["time"]["hydraulic_timestep"]
        constant = rho * gravity * time_step * inv(efficiency)

        for (a, pump) in nw_ref[:pump]
            if haskey(pump, "energy_price")
                energy_price = pump["energy_price"]
                pump_curve = ref(wm, n, :pump, a)["pump_curve"]
                curve_fun = _get_function_from_pump_curve(pump_curve)

                q, x_pump = [var(wm, n)[:q][a], var(wm, n)[:x_pump][a]]
                q_lb, q_ub = [JuMP.lower_bound(q), JuMP.upper_bound(q)]

                # Generate a set of uniform flow and cubic function breakpoints.
                breakpoints = range(q_lb, stop=q_ub, length=wm.ext[:num_breakpoints])
                f = _get_pump_gain_values(collect(breakpoints), curve_fun)

                # Add the cost corresponding to the current pump's operation.
                K = 1:wm.ext[:num_breakpoints]
                cost = constant*energy_price*sum(f[k] * lambda[a, k] for k in K)
                JuMP.add_to_expression!(objective, cost)
            else
                Memento.error(_LOGGER, "No cost given for pump \"$(pump["name"])\"")
            end
        end
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end
