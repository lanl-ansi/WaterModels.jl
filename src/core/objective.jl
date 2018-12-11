function objective_minimize_gamma(wm::GenericWaterModel)
    arcs_from = collect(ids(wm, :pipes))
    return @objective(wm.model, Min, sum(sum(wm.var[:nw][n][:gamma][a] for a in arcs_from) for n in nws(wm)))
end

function objective_minimize_cost(wm::GenericWaterModel)
    cost_function = get_diameter_cost_function(wm)
    @constraint(wm.model, wm.var[:nw][wm.cnw][:objective] >= cost_function)
    return @objective(wm.model, Min, wm.var[:nw][wm.cnw][:objective])
end

function get_diameter_cost_function(wm::GenericWaterModel)
    cost_function = zero(AffExpr)

    for n in nws(wm)
        for (a, pipe) in wm.ref[:nw][n][:ne_pipe]
            length = pipe["length"]
            diameter_vars = wm.var[:nw][n][:psi][a]
            costs = [d["costPerUnitLength"] * length for d in pipe["diameters"]]
            cost_function += AffExpr(diameter_vars[:], costs, 0.0)
        end
    end

    return cost_function
end

function objective_minimize_resistance_cost(wm::GenericWaterModel)
    objective = zero(AffExpr)
    max_cost = 0.0

    for n in nws(wm)
        C = calc_resistance_costs(wm, n)

        for (a, connection) in wm.ref[:nw][n][:connection]
            L_a = connection["length"]

            for r in 1:length(C[a])
                objective += (L_a * C[a][r]) * wm.var[:nw][n][:xr][a][r]
            end

            max_cost += L_a * C[a][1]
        end
    end

    return @objective(wm.model, Min, objective)
end

function objective_cvx_hw(wm::GenericWaterModel)
    objective = zero(AffExpr)

    for n in nws(wm)
        # TODO: Address the efficiency here, too.
        R = calc_resistances_hw(wm, n)

        for (a, connection) in wm.ref[:nw][n][:connection]
            q_p = wm.var[:nw][n][:qp][a][1]
            q_n = wm.var[:nw][n][:qn][a][1]
            L = connection["length"]
            coeff = R[a][1] * L * 0.350631
            term = @variable(wm.model, lowerbound = 0.0, category = :Cont, start = 1.0e-6)
            @NLconstraint(wm.model, term >= coeff * (q_p * (q_p^2)^0.926 + q_n * (q_n^2)^0.926))
            objective += term
        end
    end

    for n in nws(wm)
        for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
            for (a, connection) in wm.ref[:nw][n][:connection]
                q_p = wm.var[:nw][n][:qp][a][1]
                q_n = wm.var[:nw][n][:qn][a][1]
                objective -= reservoir["head"] * (q_p - q_n)
            end
        end
    end

    return @objective(wm.model, Min, objective)
end

function objective_maximize_variable(wm::GenericWaterModel, variable::JuMP.Variable)
    return @objective(wm.model, Max, variable)
end

function objective_minimize_variable(wm::GenericWaterModel, variable::JuMP.Variable)
    return @objective(wm.model, Min, variable)
end

function objective_dummy(wm::GenericWaterModel)
    return @objective(wm.model, Min, 0.0)
end
