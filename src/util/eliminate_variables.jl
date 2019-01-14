function eliminate_variables(wm::GenericWaterModel, nlp::GenericWaterModel,
                             solver::MathProgBase.AbstractMathProgSolver,
                             n::Int = wm.cnw)
    # Set common parameters.
    setsolver(nlp.model, solver)
    connection_ids = collect(ids(wm, n, :connection))
    resistance_indices = Dict{Int, Int}(a => 1 for a in connection_ids)
    num_resistances = num_resistances_eliminated = 0

    # Resistance indices.
    for (a, connection) in wm.ref[:nw][n][:connection]
        resistance_indices[a] = length(wm.ref[:nw][n][:resistance][a])
        num_resistances += length(wm.ref[:nw][n][:resistance][a])
    end

    for (a, connection) in wm.ref[:nw][n][:connection]
        for resistance_index in 1:length(wm.ref[:nw][n][:resistance][a])
            resistance_indices[a] = resistance_index
            q, h = get_cvx_solution(wm, resistance_indices, solver)
            qlb, qub, hlb, hub = check_solution_bounds(wm, q, h, resistance_indices, n)
            solution_is_feasible = all([all(values(hlb)), all(values(hub))])

            # TODO: Why does the below entail bad cuts on some instances (e.g., Hanoi)?
            #solution_is_feasible = all([all(values(qlb)), all(values(qub)),
            #                            all(values(hlb)), all(values(hub))])

            if !solution_is_feasible
                setupperbound(wm.var[:nw][n][:xr][a][resistance_index], 0)
                setupperbound.(wm.var[:nw][n][:qn][a][:, resistance_index], 0.0)
                setupperbound.(wm.var[:nw][n][:qp][a][:, resistance_index], 0.0)
                setupperbound.(wm.var[:nw][n][:xsn][a][:, resistance_index], 0)
                setupperbound.(wm.var[:nw][n][:xsp][a][:, resistance_index], 0)
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

    println("$(num_resistances_eliminated) of $(num_resistances) resistances eliminated.")
end
