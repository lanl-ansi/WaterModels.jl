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
    objective_des(wm::AbstractWaterModel)

Sets the objective function for network design (des) problem specifications. By default, the
cost of selecting the discrete pipe resistances over all design pipes is minimized.
"""
function objective_des(wm::AbstractWaterModel)
    objective = JuMP.AffExpr(0.0)

    for (n, nw_ref) in nws(wm)
        for (a, des_pipe) in ref(wm, n, :des_pipe)
            costs = des_pipe["length"] .* ref(wm, n, :resistance_cost, a)
            JuMP.add_to_expression!(objective, sum(costs .* var(wm, n, :z_des_pipe, a)))
        end
    end

    JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end
