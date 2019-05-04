function get_resistance_cost_expression(wm::GenericWaterModel)
    cost = zero(AffExpr)

    for n in nws(wm)
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
    return @objective(wm.model, Min, objective)
end

function objective_cnlp_hw(wm::GenericWaterModel)
    # Register the integrated head loss JuMP function.
    function_head_loss_integrated_hw(wm)

    # Initialize the objective.
    objective = zero(AffExpr)

    connection_ids = collect(ids(wm, :connection))

    for n in nws(wm)
        terms = @variable(wm.model, [a in connection_ids], start = 0.0,
                          lowerbound = 0.0, category = :Cont,
                          basename = "term_$(n)")

        for (a, connection) in wm.ref[:nw][n][:connection]
            q_n = wm.var[:nw][n][:qn][a][1]
            q_p = wm.var[:nw][n][:qp][a][1]

            L = connection["length"]
            coeff = wm.ref[:nw][n][:resistance][a][1] * L

            con = @NLconstraint(wm.model, coeff * (head_loss_integrated_hw(q_p) +
                                head_loss_integrated_hw(q_n)) <= terms[a])
            objective += terms[a]
        end
    end

    for n in nws(wm)
        connections = wm.ref[:nw][n][:connection]

        for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
            for (a, connection) in filter(a -> i == parse(Int, a.second["node1"]), connections)
                q_n = wm.var[:nw][n][:qn][a][1]
                q_p = wm.var[:nw][n][:qp][a][1]
                objective -= reservoir["head"] * (q_p - q_n)
            end

            for (a, connection) in filter(a -> i == parse(Int, a.second["node2"]), connections)
                q_n = wm.var[:nw][n][:qn][a][1]
                q_p = wm.var[:nw][n][:qp][a][1]
                objective -= reservoir["head"] * (q_n - q_p)
            end
        end
    end

    return @objective(wm.model, Min, objective)
end
