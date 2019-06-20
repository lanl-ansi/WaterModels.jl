######################################################################
# This file defines commonly-used objectives for water systems models.
######################################################################

function get_resistance_cost_expression(wm::GenericWaterModel, n::Int=wm.cnw)
    expr = JuMP.AffExpr(0.0)

    for (a, link) in wm.ref[:nw][n][:links_ne]
        x_res = wm.var[:nw][n][:x_res][a]
        costs = link["length"] .* wm.ref[:nw][n][:resistance_cost][a]
        num_resistances = length(x_res)
        JuMP.add_to_expression!(expr, sum(costs[r] * x_res[r] for r in 1:num_resistances))
    end

    return expr
end

function objective_ne(wm::GenericWaterModel, n::Int=wm.cnw)
    objective = get_resistance_cost_expression(wm, n)
    return JuMP.@objective(wm.model, MOI.MIN_SENSE, objective)
end
