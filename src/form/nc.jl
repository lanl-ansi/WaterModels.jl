# Define NC (nonconvex nonlinear programming) implementations of water distribution
# feasibility and optimization problem specifications. Note that here, ``nonconvex
# nonlinear'' describes the treatment of head loss and head gain constraints, which are
# nonlinear equalities. NC formulations also include nonlinearites that would be surmised in
# a full problem formulation (e.g., in the optimal water flow problem, a cubic function of
# flow is used in this formulation to define the consumption of power by active pumps).


"Adds head loss constraint for a pipe in the `NC` formulation."
function constraint_pipe_head_loss(wm::NCWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, alpha::Float64, L::Float64, r::Float64)
    # Gather flow and head variables included in head loss constraints.
    q, h_i, h_j = var(wm, n, :q_pipe, a), var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # Add nonconvex constraint for the head loss relationship.
    c = JuMP.@NLconstraint(wm.model, r * head_loss(q) == inv(L) * (h_i - h_j))

    # Append the :pipe_head_loss constraint array.
    append!(con(wm, n, :pipe_head_loss)[a], [c])
end


"Adds head loss constraint for a design pipe in the `NC` formulation."
function constraint_on_off_pipe_head_loss_des(wm::NCWaterModel, n::Int, a::Int, exponent::Float64, node_fr::Int, node_to::Int, L::Float64, resistances)
    # Gather common flow and head variables, as well as design indices.
    q, R = var(wm, n, :q_des_pipe, a), 1:length(resistances)
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # Add the nonconvex, design-expanded head loss constraint.
    lhs = JuMP.@NLexpression(wm.model, sum(resistances[r] * head_loss(q[r]) for r in R))
    c = JuMP.@NLconstraint(wm.model, lhs == inv(L) * (h_i - h_j))

    # Append the :on_off_pipe_head_loss_des constraint array.
    append!(con(wm, n, :on_off_pipe_head_loss_des)[a], [c])
end


"Adds head gain constraints for pumps in `NC` formulations."
function constraint_on_off_pump_head_gain(wm::NCWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, pc::Array{Float64}, q_min_active::Float64)
    # Gather pump flow, head gain, and status variables.
    q, g, z = var(wm, n, :q_pump, a), var(wm, n, :g_pump, a), var(wm, n, :z_pump, a)

    # Add constraint equating head gain with respect to the pump curve.
    c = JuMP.@constraint(wm.model, pc[1]*q^2 + pc[2]*q + pc[3]*z == g)

    # Append the :on_off_pump_head_gain constraint array.
    append!(con(wm, n, :on_off_pump_head_gain)[a], [c])
end


"Defines the objective for the owf problem is `NC` formulations."
function objective_owf(wm::NCWaterModel)
    objective = zero(JuMP.QuadExpr)

    for (n, nw_ref) in nws(wm)
        efficiency = 0.85 # TODO: How can the efficiency curve be used?
        rho = 1000.0 # Water density (kilogram per cubic meter).
        gravity = 9.80665 # Gravitational acceleration (meter per second squared).
        coeff = rho * gravity * ref(wm, n, :time_step) * inv(efficiency)

        for (a, pump) in nw_ref[:pump]
            if haskey(pump, "energy_price")
                # Get pump flow and head gain variables.
                q, g = var(wm, n, :q_pump, a), var(wm, n, :g_pump, a)

                # Constrain cost_var and append to the objective expression.
                JuMP.add_to_expression!(objective, coeff * pump["energy_price"], q, g)
            else
                Memento.error(_LOGGER, "No cost given for pump \"$(pump["name"])\"")
            end
        end
    end

    # Minimize the cost (in units of currency) required to operate pumps.
    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end
