######################################################################
# This file defines commonly-used objectives for water systems models.
######################################################################

function get_resistance_cost_expression(wm::GenericWaterModel, n::Int=wm.cnw)
    expr = JuMP.AffExpr(0.0)

    for (a, link) in ref(wm, n, :links_ne)
        costs = link["length"] .* ref(wm, n, :resistance_cost, a)
        JuMP.add_to_expression!(expr, sum(costs .* var(wm, n, :x_res, a)))
    end

    return expr
end

function objective_ne(wm::GenericWaterModel, n::Int=wm.cnw)
    objective = get_resistance_cost_expression(wm, n)
    return JuMP.@objective(wm.model, MOI.MIN_SENSE, objective)
end
