function find_initial_solution(wm::GenericWaterModel, params::Dict{String, Any}, nlp_solver::MathProgBase.AbstractMathProgSolver)
    n = wm.cnw
    num_rounds = 1
    connection_ids = collect(ids(wm, wm.cnw, :connection))
    R_id = Dict{Int, Int}(a => 1 for a in connection_ids)
    R_id_best = Dict{Int, Int}(a => 1 for a in connection_ids)
    params["obj_best"] = 0.0 # Initialize the best objective.

    # Solutions with minimum resistance should be automatically feasible.
    for (a, connection) in wm.ref[:nw][n][:connection]
        L_a = connection["length"]
        R_id_best[a] = length(wm.ref[:nw][n][:resistance][a])
        params["obj_best"] += L_a * wm.ref[:nw][n][:resistance_cost][a][end]
    end

    while num_rounds <= 50
        objective = 0.0

        for (a, connection) in wm.ref[:nw][n][:connection]
            L_a = connection["length"]
            R_id[a] = rand(1:length(wm.ref[:nw][n][:resistance][a]))
            objective += L_a * wm.ref[:nw][n][:resistance_cost][a][R_id[a]]
        end

        while objective > 0.85 * params["obj_best"]
            argmax_cost_reduction = connection_ids[1]
            max_cost_reduction = 0.0

            for (a, connection) in wm.ref[:nw][n][:connection]
                L_a = connection["length"]
                C_a = wm.ref[:nw][n][:resistance_cost][a]

                if R_id[a] > 1
                    cost_reduction = L_a * (C_a[R_id[a]] - C_a[R_id[a]-1])

                    if cost_reduction > max_cost_reduction
                        max_cost_reduction = cost_reduction
                        argmax_cost_reduction = a
                    end
                end
            end

            R_id[argmax_cost_reduction] -= 1
            objective -= max_cost_reduction
        end

        repaired, R_id = repair_solution(wm, R_id, params, nlp_solver)

        if repaired && objective < params["obj_best"]
            params["obj_best"] = objective
            R_id_best = deepcopy(R_id)
        end

        num_rounds += 1
    end

    return R_id_best
end
