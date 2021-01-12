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
            z_des_pipe_term = des_pipe["cost"] * var(wm, n, :z_des_pipe, a)
            JuMP.add_to_expression!(objective, z_des_pipe_term)
        end
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end


"""
    objective_owf(wm::AbstractWaterModel)

Sets the objective function for optimal water flow (owf) problem specifications.
"""
function objective_owf(wm::AbstractWaterModel)
    objective = JuMP.AffExpr(0.0)

    for (n, nw_ref) in nws(wm)
        for (a, pump) in ref(wm, n, :pump)
            @assert haskey(pump, "energy_price") # Ensure a price exists.
            coeff = _DENSITY * _GRAVITY * ref(wm, n, :time_step) * pump["energy_price"]
            JuMP.add_to_expression!(objective, coeff * var(wm, n, :Ps_pump, a))
        end
    end

    # Minimize the cost (in units of currency) required to operate pumps.
    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end