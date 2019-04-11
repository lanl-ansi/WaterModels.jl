#function objective_minimize_gamma(wm::GenericWaterModel)
#    arcs_from = collect(ids(wm, :pipes))
#    return @objective(wm.model, Min, sum(sum(wm.var[:nw][n][:gamma][a] for a in arcs_from) for n in nws(wm)))
#end
#
#function objective_minimize_cost(wm::GenericWaterModel)
#    cost_function = get_diameter_cost_function(wm)
#    @constraint(wm.model, wm.var[:nw][wm.cnw][:objective] >= cost_function)
#    return @objective(wm.model, Min, wm.var[:nw][wm.cnw][:objective])
#end

#function get_diameter_cost_function(wm::GenericWaterModel)
#    cost_function = zero(AffExpr)
#
#    for n in nws(wm)
#        for (a, pipe) in wm.ref[:nw][n][:ne_pipe]
#            length = pipe["length"]
#            diameter_vars = wm.var[:nw][n][:psi][a]
#            costs = [d["costPerUnitLength"] * length for d in pipe["diameters"]]
#            cost_function += AffExpr(diameter_vars[:], costs, 0.0)
#        end
#    end
#
#    return cost_function
#end

function get_resistance_cost_expression(wm::GenericWaterModel)
    cost = JuMP.@expression(wm.model, 0.0)

    for n in nw_ids(wm)
        for (a, connection) in wm.ref[:nw][n][:connection]
            xr = wm.var[:nw][n][:xr][a]
            C_a = connection["length"] * wm.ref[:nw][n][:resistance_cost][a]
            cost += sum(C_a[r] * xr[r] for r in 1:length(xr))
        end
    end

    return cost
end

function objective_minimize_resistance_cost(wm::GenericWaterModel)
    objective = get_resistance_cost_expression(wm)
    JuMP.@objective(wm.model, MOI.MIN_SENSE, objective)
end

#function objective_cvx_hw(wm::GenericWaterModel, n::Int = wm.cnw)
#    # Register the integrated head loss JuMP function.
#    function_head_loss_integrated_hw(wm)
#
#    # Initialize the objective.
#    objective = zero(AffExpr)
#
#    connection_ids = collect(ids(wm, :connection))
#
#    for (a, connection) in wm.ref[:nw][n][:connection]
#        q_n = wm.var[:nw][n][:q_n][a][1]
#        q_p = wm.var[:nw][n][:q_p][a][1]
#
#        L = connection["length"]
#        coeff = wm.ref[:nw][n][:resistance][a][1] * L
#        con = @NLconstraint(wm.model, coeff * (head_loss_integrated_hw(q_p) + head_loss_integrated_hw(q_n)) <= terms[a])
#
#        objective += terms[a]
#    end
#
#    for n in nws(wm)
#        connections = wm.ref[:nw][n][:connection]
#
#        for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
#            for (a, connection) in filter(a -> i == parse(Int, a.second["node1"]), connections)
#                q_n = wm.var[:nw][n][:qn][a][1]
#                q_p = wm.var[:nw][n][:qp][a][1]
#                objective -= reservoir["head"] * (q_p - q_n)
#            end
#
#            for (a, connection) in filter(a -> i == parse(Int, a.second["node2"]), connections)
#                q_n = wm.var[:nw][n][:qn][a][1]
#                q_p = wm.var[:nw][n][:qp][a][1]
#                objective -= reservoir["head"] * (q_n - q_p)
#            end
#        end
#    end
#
#    return @objective(wm.model, :Max, -objective)
#end

#function objective_maximize_variable(wm::GenericWaterModel, variable::JuMP.Variable)
#    return @objective(wm.model, Max, variable)
#end
#
#function objective_minimize_variable(wm::GenericWaterModel, variable::JuMP.Variable)
#    return @objective(wm.model, Min, variable)
#end
#
#function objective_dummy(wm::GenericWaterModel)
#    return @objective(wm.model, Min, 0.0)
#end
