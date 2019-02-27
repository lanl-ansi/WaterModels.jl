export heuristic_cut_callback_generator

import JuMP
import MathProgBase

function heuristic_cut_callback_generator(wm::GenericWaterModel, params::Dict{String, Any},
                                          nlp_solver::MathProgBase.AbstractMathProgSolver,
                                          n::Int = wm.cnw)
    function heuristic_cut_callback(cb::MathProgBase.MathProgCallbackData)
        num_nodes_explored = convert(Int, MathProgBase.cbgetexplorednodes(cb))
        resistances = wm.ref[:nw][n][:resistance]
        connection_ids = collect(ids(wm, n, :connection))
        resistance_indices = Dict{Int, Int}(a => 1 for a in connection_ids)

        if (num_nodes_explored > 0) && (num_nodes_explored % 500 == 0)
            for (a, connection) in wm.ref[:nw][n][:connection]
                xr_a = getvalue(wm.var[:nw][n][:xr][a])
                geq_indices = filter(r -> xr_a[r] >= 1.0 / length(xr_a), 1:length(xr_a))
                r_val, r_rel = findmin([xr_a[r] for r in geq_indices])
                resistance_indices[a] = geq_indices[r_rel]
            end
        else
            # Update resistances used throughout the network.
            for (a, connection) in wm.ref[:nw][n][:connection]
                xr_a = getvalue(wm.var[:nw][n][:xr][a])
                xr_ones = findall(r -> isapprox(xr_a[r], 1.0, atol = 0.01), 1:length(xr_a))
                xr_zeros = findall(r -> isapprox(xr_a[r], 0.0, atol = 0.01), 1:length(xr_a))
                xr_are_integers = (length(xr_ones) + length(xr_zeros)) == length(xr_a)

                dir = getvalue(wm.var[:nw][n][:dir][a])
                dir_is_integer = isapprox(dir, 1.0, atol = 0.01) ||
                                 isapprox(dir, 0.0, atol = 0.01)

                if !(xr_are_integers && dir_is_integer)
                    return
                else
                    resistance_indices[a] = xr_ones[1]
                end
            end
        end

        # TODO: Toggle this on and off and examine the effects.
        if compute_objective(wm, resistance_indices, n) >= params["obj_best"]
            return
        end

        repaired, resistance_indices = repair_solution(wm, resistance_indices,
                                                       params["max_repair_iters"],
                                                       params["obj_best"],
                                                       nlp_solver, n)

        if repaired
            # Update objective values.
            current_objective = compute_objective(wm, resistance_indices, n)
            params["obj_last"] = params["obj_curr"]
            params["obj_curr"] = current_objective
            params["obj_best"] = min(params["obj_curr"], params["obj_best"])

            q, h = get_cvx_solution(wm, resistance_indices, nlp_solver)

            # Set the integer values appropriately.
            for (a, connection) in wm.ref[:nw][n][:connection]
                setsolutionvalue(cb, wm.var[:nw][n][:dir][a], q[a] >= 0.0 ? 1 : 0)

                resistance_index = resistance_indices[a]
                setsolutionvalue(cb, wm.var[:nw][n][:xr][a][resistance_index], 1)

                for r in setdiff(1:length(resistances[a]), [resistance_index])
                    setsolutionvalue(cb, wm.var[:nw][n][:xr][a][r], 0)
                end

                x_par = wm.var[:nw][n][:x_par][a]
                lambda_n = wm.var[:nw][n][:lambda_n][a]
                lambda_p = wm.var[:nw][n][:lambda_p][a]

                for r in 1:length(resistances[a])
                    x_res = wm.var[:nw][n][:xr][a][r]
                    qn = wm.ref[:nw][n][:qn_bp][a][r]
                    qp = wm.ref[:nw][n][:qp_bp][a][r]

                    for k in 1:length(qp)-1
                        if q[a] > 0.0 && r == resistance_indices[a]
                            if q[a] >= qp[k] && q[a] <= qp[k+1]
                                setsolutionvalue(cb, x_par[k], 1)
                                lambda_val = (q[a] - qp[k]) / (qp[k+1] - qp[k])
                                setsolutionvalue(cb, lambda_p[r, k+1], lambda_val)
                                setsolutionvalue(cb, lambda_p[r, k], 1.0 - lambda_val)
                            end
                        elseif q[a] <= 0.0 && r == resistance_indices[a]
                            if -q[a] >= qn[k] && -q[a] <= qn[k+1]
                                setsolutionvalue(cb, x_par[k], 1)
                                lambda_val = (-q[a] - qn[k]) / (qn[k+1] - qn[k])
                                setsolutionvalue(cb, lambda_n[r, k+1], lambda_val)
                                setsolutionvalue(cb, lambda_n[r, k], 1.0 - lambda_val)
                            end
                        end
                    end
                end
            end

            # Register the solution via the callback.
            addsolution(cb)
        end
    end

    return heuristic_cut_callback
end
