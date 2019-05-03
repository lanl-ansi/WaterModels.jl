function eliminate_variables(wm::GenericWaterModel, nlp::GenericWaterModel,
                             solver::MathProgBase.AbstractMathProgSolver,
                             n::Int = wm.cnw)
    # Set common parameters.
    setsolver(nlp.model, solver)
    resistances = wm.ref[:nw][n][:resistance]
    connection_ids = collect(ids(wm, n, :connection))

    # Initialize resistances.
    resistance_indices = Dict{Int, Int}(a => length(resistances[a]) for a in connection_ids)
    num_resistances = sum([length(resistances[a]) for a in connection_ids])
    num_resistances_eliminated = 0

    for (a, connection) in wm.ref[:nw][n][:connection]
        for resistance_index in 1:length(resistances[a])
            resistance_indices[a] = resistance_index
            q, h = get_cnlp_solution(wm, resistance_indices, solver)
            qlb, qub, hlb, hub = check_solution_bounds(wm, q, h, resistance_indices, n)

            # TODO: Why does the below entail bad cuts on some instances (e.g., Hanoi)?
            #solution_is_feasible = all([all(values(qlb)), all(values(qub)),
            #                            all(values(hlb)), all(values(hub))])

            solution_is_feasible = all([all(values(hlb)), all(values(hub))])

            if !solution_is_feasible
                setupperbound(wm.var[:nw][n][:xr][a][resistance_index], 0)
                setupperbound(wm.var[:nw][n][:qn][a][resistance_index], 0.0)
                setupperbound(wm.var[:nw][n][:qp][a][resistance_index], 0.0)
                setupperbound(nlp.var[:nw][n][:xr][a][resistance_index], 0)
                setupperbound(nlp.var[:nw][n][:qn][a][resistance_index], 0.0)
                setupperbound(nlp.var[:nw][n][:qp][a][resistance_index], 0.0)
                num_resistances_eliminated += 1
            else
                break
            end
        end

        resistance_indices[a] = length(wm.ref[:nw][n][:resistance][a])
    end
end
