# Define NCNLP (nonconvex nonlinear programming) implementations of water
# distribution feasibility and optimization problem specifications. Note that
# here, ``nonconvex nonlinear'' describes the treatment of head loss and head
# gain constraints, which are nonlinear equalities. NCNLP formulations also
# include nonlinearites that would be surmised in a full problem formulation
# (e.g., in the optimal power flow problem, a cubic function of flow is used in
# this formulation to define the consumption of power by active pumps).

"Creates flow variables for `NCNLP` formulations (`q`)."
function variable_flow(wm::AbstractNCNLPModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    variable_flow_common(wm, nw=nw, bounded=bounded, report=report)
end

"Creates network design flow variables for `NCNLP` formulations (`q_des`, `x_des`)."
function variable_flow_des(wm::AbstractNCNLPModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    variable_flow_des_common(wm, nw=nw, bounded=bounded, report=report)
end

"Adds head loss constraints for check valves in `NCNLP` formulations."
function constraint_check_valve_head_loss(wm::AbstractNCNLPModel, n::Int, a::Int, node_fr::Int, node_to::Int, L::Float64, r::Float64) 
    # Gather common variables and data.
    q, x_cv = [var(wm, n, :q, a), var(wm, n, :x_cv, a)]
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]
    dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)

    # There are two possibilities, here: (i) x_cv = 1, in which the check valve
    # is open, and the head loss relationship should be met with equality. In
    # this case, the constraints below reduce to lhs = 0.0, as would be
    # expected. Otherwise, (ii) x_cv = 0, and the head difference must be
    # negative and decoupled from the traditional head loss relationship.
    lhs = JuMP.@NLexpression(wm.model, inv(L) * (h_i - h_j) - r * head_loss(q))
    c_1 = JuMP.@NLconstraint(wm.model, lhs <= 0.0)
    c_2 = JuMP.@NLconstraint(wm.model, lhs >= inv(L) * (1.0 - x_cv) * dh_lb)

    # Append the :head_loss constraint array.
    append!(con(wm, n, :head_loss)[a], [c_1, c_2])
end

"Adds head loss constraints for shutoff valves in `NCNLP` formulations."
function constraint_sv_head_loss(wm::AbstractNCNLPModel, n::Int, a::Int, node_fr::Int, node_to::Int, L::Float64, r::Float64) 
    # Gather common variables and data.
    q, x_sv = var(wm, n, :q, a), var(wm, n, :x_sv, a)
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)
    dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    dh_ub = JuMP.upper_bound(h_i) - JuMP.lower_bound(h_j)

    # There are two possibilities, here: (i) x_sv = 1, in which the shutoff
    # valve is open, and the head loss relationship should be met with
    # equality. In this case, the constraints below reduce to lhs = 0.0, as
    # would be expected. Otherwise, (ii) x_sv = 0, and the head difference must
    # be decoupled from the traditional head loss relationship.
    lhs = JuMP.@NLexpression(wm.model, inv(L) * (h_i - h_j) - r * head_loss(q))
    c_1 = JuMP.@NLconstraint(wm.model, lhs <= inv(L) * (1.0 - x_sv) * dh_ub)
    c_2 = JuMP.@NLconstraint(wm.model, lhs >= inv(L) * (1.0 - x_sv) * dh_lb)

    # Append the :head_loss constraint array.
    append!(con(wm, n, :head_loss)[a], [c_1, c_2])
end

"Adds head loss constraints for pipes (without check valves) in `NCNLP` formulations."
function constraint_head_loss_pipe(wm::AbstractNCNLPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    # Gather flow and head variables included in head loss constraints.
    q = var(wm, n, :q, a)
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]

    # Add nonconvex constraint for the head loss relationship.
    c = JuMP.@NLconstraint(wm.model, r * head_loss(q) == inv(L) * (h_i - h_j))

    # Append the :head_loss constraint array.
    append!(con(wm, n, :head_loss)[a], [c])
end

"Add head loss constraints for design pipes (without check valves) in `NCNLP` formulations."
function constraint_head_loss_pipe_des(wm::AbstractNCNLPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, res)
    # Gather common flow and head variables, as well as design indices.
    q, R = [var(wm, n, :q_des, a), 1:length(res)]
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]

    # Add the nonconvex, design-expanded head loss constraint.
    lhs = JuMP.@NLexpression(wm.model, sum(res[r] * head_loss(q[r]) for r in R))
    c = JuMP.@NLconstraint(wm.model, lhs == inv(L) * (h_i - h_j))

    # Append the :head_loss constraint array.
    append!(con(wm, n, :head_loss)[a], [c])
end

"Adds head gain constraints for pumps in `NCNLP` formulations."
function constraint_pump_head_gain(wm::AbstractNCNLPModel, n::Int, a::Int, node_fr::Int, node_to::Int, pc::Array{Float64})
    # Gather common flow and pump variables.
    x_pump = var(wm, n, :x_pump, a)
    q, g_var = [var(wm, n, :q, a), var(wm, n, :g, a)]

    # Add constraint equating head gain with respect to the pump curve.
    g_expr = pc[1]*q^2 + pc[2]*q + pc[3]*x_pump
    c = JuMP.@constraint(wm.model, g_expr == g_var)

    # Append the :head_gain constraint array.
    append!(con(wm, n, :head_gain)[a], [c])
end

"Defines the objective for the owf problem is `NCNLP` formulations."
function objective_owf(wm::AbstractNCNLPModel) 
    objective = JuMP.AffExpr(0.0)
    time_step = wm.ref[:option]["time"]["hydraulic_timestep"]

    for (n, nw_ref) in nws(wm)
        efficiency = 0.85 # TODO: How can the efficiency curve be used?
        rho = 1000.0 # Water density (kilogram per cubic meter).
        gravity = 9.80665 # Gravitational acceleration (meter per second squared).
        coeff = rho * gravity * time_step * inv(efficiency)

        for (a, pump) in nw_ref[:pump]
            if haskey(pump, "energy_price")
                # Get price data and power-related variables.
                price = pump["energy_price"]
                q, g = [var(wm, n)[:q][a], var(wm, n)[:g][a]]

                # Constrain cost_var and append to the objective expression.
                cost_var = JuMP.@variable(wm.model, lower_bound=0.0)
                cost_expr = JuMP.@NLexpression(wm.model, coeff * price * g * q)
                c = JuMP.@NLconstraint(wm.model, cost_expr <= cost_var)
                JuMP.add_to_expression!(objective, cost_var)
            else
                Memento.error(_LOGGER, "No cost given for pump \"$(pump["name"])\"")
            end
        end
    end

    # Minimize the cost (in units of currency) required to operate pumps.
    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end
