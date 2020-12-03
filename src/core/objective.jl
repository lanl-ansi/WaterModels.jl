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

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end


function objective_owf_best_efficiency(wm::AbstractWaterModel)
    objective = JuMP.AffExpr(0.0)

    for (n, nw_ref) in nws(wm)
        for (a, pump) in nw_ref[:pump]
            if haskey(pump, "energy_price")
                # Get pump best efficiency data required for construction.
                flow = _calc_pump_best_efficiency_flow(pump)
                power = _calc_pump_best_efficiency_power(pump)

                # Get constant data associated with the cost function.
                constant = pump["energy_price"] * ref(wm, n, :time_step) * power

                # Get flow-related variables and data.
                q, z = var(wm, n, :q_pump, a), var(wm, n, :z_pump, a)

                # Add the cost corresponding to the current pump's operation.
                cost = constant * (inv(3.0) * q * inv(flow) + z * 2.0 * inv(3.0))
                JuMP.add_to_expression!(objective, cost)
            else
                Memento.error(_LOGGER, "No cost given for pump \"$(pump["name"])\"")
            end
        end
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end


"""
    objective_owf(wm::AbstractWaterModel)

Sets the objective function for optimal water flow (owf) problem specifications.
"""
function objective_owf(wm::AbstractWaterModel)
    if get(wm.ext, :use_best_efficiency_form, false)
        return objective_owf_best_efficiency(wm)
    else
        return objective_owf_default(wm)
    end
end
