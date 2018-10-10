function objective_minimize_gamma(wm::GenericWaterModel)
    arcs_from = collect(ids(wm, :pipes))
    return @objective(wm.model, Min, sum(sum(wm.var[:nw][n][:gamma][a] for a in arcs_from) for n in nws(wm)))
end

function objective_minimize_cost(wm::GenericWaterModel)
    cost_function = get_diameter_cost_function(wm)
    @constraint(wm.model, wm.var[:nw][wm.cnw][:objective] == cost_function)
    return @objective(wm.model, Min, wm.var[:nw][wm.cnw][:objective])
    #return @NLobjective(wm.model, Min, wm.var[:nw][wm.cnw][:objective])
end

function get_diameter_cost_function(wm::GenericWaterModel)
    cost_function = @expression(wm.model, 0.0)

    for n in nws(wm)
        for (a, pipe) in wm.ref[:nw][n][:ne_pipe]
            length = pipe["length"]
            diameter_vars = wm.var[:nw][n][:psi][a]
            diameters = [key[1] for key in keys(diameter_vars)]
            costs = [d["costPerUnitLength"] * length for d in pipe["diameters"]] * 1.0e-6
            cost_function += AffExpr(diameter_vars[:], costs, 0.0)
        end
    end
end

function objective_maximize_variable(wm::GenericWaterModel, variable::JuMP.Variable)
    return @objective(wm.model, Max, variable)
end

function objective_maximize_variable(wm::GenericWaterModel, variable::JuMP.Variable)
    return @objective(wm.model, Min, variable)
end

function objective_dummy(wm::GenericWaterModel)
    return @NLobjective(wm.model, Min, 0.0)
    #return @objective(wm.model, Min, 0.0)
end
