# Define NCNLP (non-convex nonlinear programming) implementations of water distribution models.

function constraint_head_loss_check_valve(wm::AbstractNCNLPModel, n::Int, a::Int, node_fr::Int, node_to::Int, L::Float64, r::Float64) 
    q = var(wm, n, :q, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)
    x_cv = var(wm, n, :x_cv, a)

    ## TODO: Possible formulation below (expand out...)?
    #c = JuMP.@NLconstraint(wm.model, x_cv * r * head_loss(q) == inv(L) * x_cv * (h_i - h_j))
    #c = JuMP.@NLconstraint(wm.model, (1 - x_cv) * h_j >= (1 - x_cv) * h_i)

    lhs = JuMP.@NLexpression(wm.model, r * head_loss(q) - inv(L) * (h_i - h_j))
    c_1 = JuMP.@NLconstraint(wm.model, lhs <= 1.0e6 * (1.0 - x_cv))
    c_2 = JuMP.@NLconstraint(wm.model, lhs >= -1.0e6 * (1.0 - x_cv))
    append!(con(wm, n, :head_loss)[a], [c_1, c_2])
end

function constraint_head_loss_pipe(wm::AbstractNCNLPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, r::Float64) 
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)
    q = var(wm, n, :q, a)

    c = JuMP.@NLconstraint(wm.model, L*r * head_loss(q) == (h_i - h_j))
    append!(con(wm, n, :head_loss)[a], [c])
end

function constraint_head_loss_pipe_ne(wm::AbstractNCNLPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, pipe_resistances)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)
    q_ne = var(wm, n, :q_ne, a)

    lhs = JuMP.@NLexpression(wm.model, sum(L*r * head_loss(q_ne[r_id]) for
        (r_id, r) in enumerate(pipe_resistances)))
    c = JuMP.@NLconstraint(wm.model, lhs == h_i - h_j)
    append!(con(wm, n, :head_loss)[a], [c])
end

"Pump head gain constraint when the pump status is ambiguous."
function constraint_head_gain_pump(wm::AbstractNCNLPModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Gather common variables.
    q = var(wm, n, :q, a)
    g = var(wm, n, :g, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)
    x_pump = var(wm, n, :x_pump, a)
    q_ub = JuMP.has_upper_bound(q) ? JuMP.upper_bound(q) : 10.0

    # Define the head gain caused by the pump.
    c_1 = JuMP.@NLconstraint(wm.model, curve_fun[1]*q^2 + curve_fun[2]*q + curve_fun[3]*x_pump == g)

    # If the pump is off, decouple the head difference relationship.
    c_2 = JuMP.@constraint(wm.model, (h_j - h_i) - g <= 1.0e6 * (1.0 - x_pump))
    c_3 = JuMP.@constraint(wm.model, (h_j - h_i) - g >= -1.0e6 * (1.0 - x_pump))

    # If the pump is off, the flow along the pump must be zero.
    c_4 = JuMP.@constraint(wm.model, q <= q_ub * x_pump)
    c_5 = JuMP.@constraint(wm.model, q >= 1.0e-6 * x_pump)

    # Append the constraint array.
    append!(con(wm, n, :head_gain)[a], [c_1, c_2, c_3, c_4, c_5])
end

"Pump head gain constraint when the pump is forced to be on."
function constraint_head_gain_pump_on(wm::AbstractNCNLPModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Gather common variables.
    q = var(wm, n, :q, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)

    # Define the head difference relationship when the pump is on (h_j >= h_i).
    c = JuMP.@NLconstraint(wm.model, curve_fun[1]*q^2 + curve_fun[2]*q + curve_fun[3] == (h_j - h_i))
    con(wm, n, :head_gain)[a] = [c]
end

function objective_owf(wm::AbstractNCNLPModel) 
    expr = JuMP.AffExpr(0.0)

    for (n, nw_ref) in nws(wm)
        pump_ids = ids(wm, n, :pump)
        costs = JuMP.@variable(wm.model, [a in pump_ids], base_name="costs[$(n)]",
                               start=get_start(ref(wm, n, :pump), a, "costs", 0.0))

        efficiency = 0.85 # TODO: Change this after discussion. 0.85 follows Fooladivanda.
        rho = 1000.0 # Water density.
        g = 9.80665 # Gravitational acceleration.

        time_step = nw_ref[:option]["time"]["hydraulic_timestep"]
        constant = rho * g * time_step * inv(efficiency)

        for (a, pump) in nw_ref[:pump]
            if haskey(pump, "energy_price")
                g = var(wm, n)[:g][a]
                q = var(wm, n)[:q][a]

                energy_price = pump["energy_price"]
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
