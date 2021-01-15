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
            coeff = ref(wm, n, :time_step) * pump["energy_price"] # * _DENSITY * _GRAVITY
            JuMP.add_to_expression!(objective, coeff * var(wm, n, :Ps_pump, a))
        end
    end

    # Normalize the objective to be more numerically well-behaved.
    pos_coeff = filter(x -> x > 0.0, collect(objective.terms.vals))
    minimum_scalar = length(pos_coeff) > 0 ? minimum(pos_coeff) : 1.0
    objective_scaled = (1.0 / minimum_scalar) * objective

    # Minimize the (numerically scaled) cost required to operate pumps.
    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective_scaled)
end