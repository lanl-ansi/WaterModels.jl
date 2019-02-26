function bound_tightening(wm::GenericWaterModel,
                          nlp::GenericWaterModel,
                          best_objective_value::Float64,
                          num_iterations::Int,
                          solver::MathProgBase.AbstractMathProgSolver,
                          n::Int = wm.cnw)
    # Add a constraint to restrict the cost.
    cost = get_resistance_cost_expression(nlp)
    @constraint(nlp.model, cost <= best_objective_value)
    setsolver(nlp.model, solver)
    connection_ids = collect(ids(wm, n, :connection))

    for num_iteration in 1:num_iterations
        # Tighten direction bounds.
        for (a, connection) in wm.ref[:nw][n][:connection]
            if sum(getupperbound.(nlp.var[:nw][n][:qn][a])) < 1.0e-9
                setlowerbound(wm.var[:nw][n][:dir][a], 1)
                setupperbound(wm.var[:nw][n][:dhn][a], 0.0)
                setlowerbound(nlp.var[:nw][n][:dir][a], 1)
                setupperbound(nlp.var[:nw][n][:dhn][a], 0.0)
            end

            if sum(getupperbound.(nlp.var[:nw][n][:qp][a])) < 1.0e-9
                setupperbound(wm.var[:nw][n][:dir][a], 0)
                setupperbound(wm.var[:nw][n][:dhp][a], 0.0)
                setupperbound(nlp.var[:nw][n][:dir][a], 0)
                setupperbound(nlp.var[:nw][n][:dhp][a], 0.0)
            end
        end

        # Update resistances used throughout the network.
        for (a, connection) in wm.ref[:nw][n][:connection]
            @objective(nlp.model, Max, sum(nlp.var[:nw][n][:qn][a]))
            status = JuMP.solve(nlp.model, relaxation = true)
            qn_objective_val = getobjectivevalue(nlp.model)

            for r in 1:length(wm.ref[:nw][n][:resistance][a])
                q_n_wm = wm.var[:nw][n][:qn][a][r]
                q_n_nlp = nlp.var[:nw][n][:qn][a][r]
                setupperbound(q_n_wm, min(getupperbound(q_n_nlp), qn_objective_val + 1.0e-6))
                setupperbound(q_n_nlp, min(getupperbound(q_n_nlp), qn_objective_val + 1.0e-6))
            end

            @objective(nlp.model, Max, sum(nlp.var[:nw][n][:qp][a]))
            status = JuMP.solve(nlp.model, relaxation = true)
            qp_objective_val = getobjectivevalue(nlp.model)

            for r in 1:length(wm.ref[:nw][n][:resistance][a])
                q_p_wm = wm.var[:nw][n][:qp][a][r]
                q_p_nlp = nlp.var[:nw][n][:qp][a][r]
                setupperbound(q_p_wm, min(getupperbound(q_p_nlp), qp_objective_val + 1.0e-6))
                setupperbound(q_p_nlp, min(getupperbound(q_p_nlp), qp_objective_val + 1.0e-6))
            end
        end

        for (i, junction) in wm.ref[:nw][n][:junctions]
            h_wm = wm.var[:nw][n][:h][i]
            h_nlp = nlp.var[:nw][n][:h][i]

            @objective(nlp.model, Min, h_nlp)
            status = JuMP.solve(nlp.model, relaxation = true)
            setlowerbound(h_wm, max(junction["elev"], getvalue(h_nlp) - 1.0e-6))
            setlowerbound(h_nlp, max(junction["elev"], getvalue(h_nlp) - 1.0e-6))

            @objective(nlp.model, Max, h_nlp)
            status = JuMP.solve(nlp.model, relaxation = true)
            setupperbound(h_wm, max(getupperbound(h_nlp), getvalue(h_nlp) + 1.0e-6))
            setupperbound(h_nlp, max(getupperbound(h_nlp), getvalue(h_nlp) + 1.0e-6))
        end
    end
end
