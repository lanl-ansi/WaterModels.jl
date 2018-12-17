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
    resistance_indices = Dict{Int, Int}(a => 1 for a in connection_ids)

    # Resistance indices.
    for (a, connection) in wm.ref[:nw][n][:connection]
        resistance_index = length(wm.ref[:nw][n][:resistance][a])
        resistance_indices[a] = resistance_index
    end

    for (a, connection) in wm.ref[:nw][n][:connection]
        for resistance_index in 1:length(wm.ref[:nw][n][:resistance][a])
            resistance_indices[a] = resistance_index
            q, h = get_cvx_solution(wm, resistance_indices, solver)
            qlb, qub, hlb, hub = check_solution_bounds(wm, q, h, resistance_indices, n)
            solution_is_feasible = all([all(values(qlb)), all(values(qub)),
                                        all(values(hlb)), all(values(hub))])

            if !solution_is_feasible
                setlowerbound(wm.var[:nw][n][:xr][a][resistance_index], 0)
                setupperbound(wm.var[:nw][n][:xr][a][resistance_index], 0)
                setlowerbound(wm.var[:nw][n][:qp][a][resistance_index], 0.0)
                setupperbound(wm.var[:nw][n][:qp][a][resistance_index], 0.0)
                setlowerbound(wm.var[:nw][n][:qn][a][resistance_index], 0.0)
                setupperbound(wm.var[:nw][n][:qn][a][resistance_index], 0.0)

                setlowerbound(nlp.var[:nw][n][:xr][a][resistance_index], 0)
                setupperbound(nlp.var[:nw][n][:xr][a][resistance_index], 0)
                setlowerbound(nlp.var[:nw][n][:qp][a][resistance_index], 0.0)
                setupperbound(nlp.var[:nw][n][:qp][a][resistance_index], 0.0)
                setlowerbound(nlp.var[:nw][n][:qn][a][resistance_index], 0.0)
                setupperbound(nlp.var[:nw][n][:qn][a][resistance_index], 0.0)
            else
                break
            end
        end

        resistance_indices[a] = length(wm.ref[:nw][n][:resistance][a])
    end

    for num_iteration in 1:num_iterations
        # Tighten direction bounds.
        for (a, connection) in wm.ref[:nw][n][:connection]
            if sum(getupperbound.(nlp.var[:nw][n][:qn][a])) < 1.0e-6
                setlowerbound(wm.var[:nw][n][:dir][a], 1)
                setlowerbound(nlp.var[:nw][n][:dir][a], 1)
                setupperbound(nlp.var[:nw][n][:dhn][a], 0.0)
                setupperbound(wm.var[:nw][n][:dhn][a], 0.0)
            end

            if sum(getupperbound.(nlp.var[:nw][n][:qp][a])) < 1.0e-6
                setupperbound(wm.var[:nw][n][:dir][a], 0)
                setupperbound(nlp.var[:nw][n][:dir][a], 0)
                setupperbound(nlp.var[:nw][n][:dhp][a], 0.0)
                setupperbound(wm.var[:nw][n][:dhp][a], 0.0)
            end
        end

        # Update resistances used throughout the network.
        for (a, connection) in wm.ref[:nw][n][:connection]
            for r in 1:length(wm.ref[:nw][n][:resistance][a])
                q_p_wm = wm.var[:nw][n][:qp][a][r]
                q_p_nlp = nlp.var[:nw][n][:qp][a][r]
                @objective(nlp.model, Max, q_p_nlp)
                status = JuMP.solve(nlp.model, relaxation = true)
                setupperbound(q_p_wm, getvalue(q_p_nlp))
                setupperbound(q_p_nlp, getvalue(q_p_nlp))

                q_n_wm = wm.var[:nw][n][:qn][a][r]
                q_n_nlp = nlp.var[:nw][n][:qn][a][r]
                @objective(nlp.model, Max, q_n_nlp)
                status = JuMP.solve(nlp.model, relaxation = true)
                setupperbound(q_n_wm, getvalue(q_n_nlp))
                setupperbound(q_n_nlp, getvalue(q_n_nlp))
            end
        end

        for (i, junction) in wm.ref[:nw][n][:junctions]
            h_wm = wm.var[:nw][n][:h][i]
            h_nlp = nlp.var[:nw][n][:h][i]

            @objective(nlp.model, Min, h_nlp)
            status = JuMP.solve(nlp.model, relaxation = true)
            setlowerbound(h_wm, getvalue(h_nlp))
            setlowerbound(h_nlp, getvalue(h_nlp))

            @objective(nlp.model, Max, h_nlp)
            status = JuMP.solve(nlp.model, relaxation = true)
            setupperbound(h_wm, getvalue(h_nlp))
            setupperbound(h_nlp, getvalue(h_nlp))
        end
    end
end
