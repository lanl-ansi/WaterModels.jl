# Define NCNLP (nonconvex nonlinear programming) implementations of water
# distribution feasibility and optimization problem specifications.

function constraint_head_loss_check_valve(wm::AbstractNCNLPModel, n::Int, a::Int, node_fr::Int, node_to::Int, L::Float64, r::Float64) 
    # Gather common variables.
    q = var(wm, n, :q, a)
    x_cv = var(wm, n, :x_cv, a)
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]

    # Add the potentially-decoupled head loss constraints.
    lhs = JuMP.@NLexpression(wm.model, inv(L) * (h_i - h_j) - r * head_loss(q))
    c_1 = JuMP.@NLconstraint(wm.model, lhs <= 0.0)
    dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    c_2 = JuMP.@NLconstraint(wm.model, lhs >= inv(L) * (1.0 - x_cv) * dh_lb)
    append!(con(wm, n, :head_loss)[a], [c_1, c_2])
end

function constraint_head_loss_pipe(wm::AbstractNCNLPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    q = var(wm, n, :q, a)
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]
    c = JuMP.@NLconstraint(wm.model, r * head_loss(q) == inv(L) * (h_i - h_j))
    append!(con(wm, n, :head_loss)[a], [c])
end

function constraint_head_loss_pipe_ne(wm::AbstractNCNLPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, pipe_resistances)
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]
    q_ne = var(wm, n, :q_ne, a)

    lhs = JuMP.@NLexpression(wm.model, sum(r * head_loss(q_ne[r_id]) for
        (r_id, r) in enumerate(pipe_resistances)))

    c = JuMP.@NLconstraint(wm.model, lhs == inv(L) * (h_i - h_j))
    append!(con(wm, n, :head_loss)[a], [c])
end

function get_curve_cut(curve_fun::Array{Float64}, q::JuMP.VariableRef, x::JuMP.VariableRef, q_min::Float64, q_max::Float64)
    g_1 = curve_fun[1]*q_min*q_min + curve_fun[2]*q_min + curve_fun[3]
    g_2 = curve_fun[1]*q_max*q_max + curve_fun[2]*q_max + curve_fun[3]
    m = (g_2 - g_1) / (q_max - q_min)
    return m * q + (g_1 - m*q_min) * x
end

"Pump head gain constraint when the pump status is ambiguous."
function constraint_head_gain_pump(wm::AbstractNCNLPModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Gather common variables.
    x_pump = var(wm, n, :x_pump, a)
    q, g = [var(wm, n, :q, a), var(wm, n, :g, a)]
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]

    # Add constraint equating the head gain variable.
    g_expr = curve_fun[1]*q*q + curve_fun[2]*q + curve_fun[3]*x_pump
    c_1 = JuMP.@constraint(wm.model, g_expr == g)

    # If the pump is off, decouple the head difference relationship.
    dh_lb = JuMP.lower_bound(h_j) - JuMP.upper_bound(h_i)
    c_3 = JuMP.@constraint(wm.model, (h_j - h_i) - g >= dh_lb * (1.0 - x_pump))
    dh_ub = JuMP.upper_bound(h_j) - JuMP.lower_bound(h_i)
    c_4 = JuMP.@constraint(wm.model, (h_j - h_i) - g <= dh_ub * (1.0 - x_pump))

    # If the pump is off, the flow along the pump must be zero.
    c_5 = JuMP.@constraint(wm.model, q <= JuMP.upper_bound(q) * x_pump)
    c_6 = JuMP.@constraint(wm.model, q >= 6.31465679e-6 * x_pump)

    # Append the constraint array.
    append!(con(wm, n, :head_gain)[a], [c_1, c_3, c_4, c_5, c_6])
end

"Pump head gain constraint when the pump is forced to be on."
function constraint_head_gain_pump_on(wm::AbstractNCNLPModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Gather common variables.
    q = var(wm, n, :q, a)
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]

    # Define the head difference relationship when the pump is on (h_j >= h_i).
    lhs = curve_fun[1]*q*q + curve_fun[2]*q + curve_fun[3]
    c = JuMP.@constraint(wm.model, lhs == h_j - h_i)
    con(wm, n, :head_gain)[a] = [c]
end

"Objective for the nonconvex optimal water flow problem."
function objective_owf(wm::AbstractNCNLPModel) 
    expr = JuMP.AffExpr(0.0)

    for (n, nw_ref) in nws(wm)
        pump_ids = ids(wm, n, :pump)
        costs = JuMP.@variable(wm.model, [a in pump_ids], base_name="costs[$(n)]",
                               lower_bound=0.0)
                               #start=get_start(ref(wm, n, :pump), a, "costs", 0.0))

        efficiency = 0.85 # TODO: Change this after discussion. 0.85 follows Fooladivanda.
        rho = 1000.0 # Water density (kilogram per cubic meter).
        gravity = 9.80665 # Gravitational acceleration (meter per second squared).

        time_step = nw_ref[:option]["time"]["hydraulic_timestep"]
        constant = rho * gravity * time_step * inv(efficiency)

        for (a, pump) in nw_ref[:pump]
            if haskey(pump, "energy_price")
                energy_price = pump["energy_price"]
                g = var(wm, n)[:g][a]
                q = var(wm, n)[:q][a]

                pump_cost = JuMP.@NLexpression(wm.model, constant * energy_price * g * q)
                con = JuMP.@NLconstraint(wm.model, pump_cost <= costs[a])
            else
                Memento.error(_LOGGER, "No cost given for pump \"$(pump["name"])\"")
            end
        end

        expr += sum(costs)
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, expr)
end
