######################################################################
# This file defines commonly-used objectives for water systems models.
######################################################################

function get_resistance_cost_expression(wm::GenericWaterModel, n::Int=wm.cnw)
    expression = JuMP.AffExpr(0.0)

    for (a, link) in wm.ref[:nw][n][:links_ne]
        xʳᵉˢ = wm.var[:nw][n][:xʳᵉˢ][a]
        costs = link["length"] .* wm.ref[:nw][n][:resistance_cost][a]
        num_resistances = length(xʳᵉˢ)
        JuMP.add_to_expression!(expression, sum(costs[r] * xʳᵉˢ[r] for r in 1:num_resistances))
    end

    return expression
end

function objective_ne(wm::GenericWaterModel, n::Int=wm.cnw)
    objective = get_resistance_cost_expression(wm, n)
    return JuMP.@objective(wm.model, MOI.MIN_SENSE, objective)
end
