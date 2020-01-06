######################################################################
# This file defines commonly-used objectives for water systems models.
######################################################################

"""
    objective_wf(wm::AbstractWaterModel)

Sets the objective function for [Water Flow (WF)](@ref) problem specifications.
By default, only feasibility must be satisfied.
"""
function objective_wf(wm::AbstractWaterModel)
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end

"""
    objective_cwf(wm::AbstractWaterModel)

Sets the objective function for constrained water flow (cwf) problem
specifications. By default, only feasibility must be satisfied.
"""
function objective_cwf(wm::AbstractWaterModel)
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end

"""
    objective_owf(wm::AbstractWaterModel)

Sets the objective function for optimal water flow (owf) problem
specifications. By default, only feasibility must be satisfied.
"""
function objective_owf(wm::AbstractWaterModel) 
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end

"""
    objective_ne(wm::AbstractWaterModel)

Sets the objective function for network expansion (ne) problem specifications.
By default, the cost of selecting discrete network resistances is minimized.
"""
function objective_ne(wm::AbstractWaterModel)
    obj = JuMP.AffExpr(0.0)

    for n in nw_ids(wm)
        for (a, link) in ref(wm, n, :link_ne)
            costs = link["length"] .* ref(wm, n, :resistance_cost, a)
            JuMP.add_to_expression!(obj, sum(costs .* var(wm, n, :x_res, a)))
        end
    end

    JuMP.@objective(wm.model, _MOI.MIN_SENSE, obj)
end
