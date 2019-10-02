######################################################################
# This file defines commonly-used objectives for water systems models.
######################################################################

function objective_wf(wm::AbstractWaterModel)
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end

function objective_cwf(wm::AbstractWaterModel)
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end

function objective_owf(wm::AbstractWaterModel) 
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end

function get_ne_expression(wm::AbstractWaterModel)
    expr = JuMP.AffExpr(0.0)

    for n in nw_ids(wm)
        for (a, link) in ref(wm, n, :link_ne)
            costs = link["length"] .* ref(wm, n, :resistance_cost, a)
            JuMP.add_to_expression!(expr, sum(costs .* var(wm, n, :x_res, a)))
        end
    end

    return expr
end

function objective_ne(wm::AbstractWaterModel)
    objective = get_ne_expression(wm)
    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end
