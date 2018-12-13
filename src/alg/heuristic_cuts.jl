export heuristic_cut_callback_generator

import JuMP
import MathProgBase

function heuristic_cut_callback_generator(wm::GenericWaterModel, params::Dict{String, Any},
                                          nlp_solver::MathProgBase.AbstractMathProgSolver,
                                          n::Int = wm.cnw)
    resistances = wm.ref[:nw][n][:resistance]
    network = deepcopy(wm.data)
    connection_ids = collect(ids(wm, n, :connection))
    R_id = Dict{Int, Int}(a => 1 for a in connection_ids)

    function heuristic_cut_callback(cb::MathProgBase.MathProgCallbackData)
        # Set up variable arrays that will be used for cuts.
        xr_ones = Array{JuMP.Variable, 1}()
        xr_zeros = Array{JuMP.Variable, 1}()

        # Initialize the objective value.
        current_objective = 0.0

        # Update resistances used throughout the network.
        for (a, connection) in wm.ref[:nw][n][:connection]
            xr_a = getvalue(wm.var[:nw][n][:xr][a])
            xr_ones = findall(r -> isapprox(xr_a[r], 1.0, atol = 0.01), 1:length(xr_a))
            xr_zeros = findall(r -> isapprox(xr_a[r], 0.0, atol = 0.01), 1:length(xr_a))
            xr_integers = (length(xr_ones) + length(xr_zeros)) == length(xr_a)

            dir = getvalue(wm.var[:nw][n][:dir][a])
            dir_integer = isapprox(dir, 1.0, atol = 0.01) || isapprox(dir, 0.0, atol = 0.01)

            if !(xr_integers && dir_integer)
                return
            else
                R_id[a] = xr_ones[1]
                L_a = connection["length"]
                network["pipes"][string(a)]["resistance"] = resistances[a][R_id[a]]
                current_objective += L_a * wm.ref[:nw][n][:resistance_cost][a][R_id[a]]
            end
        end

        repaired, R_id = repair_solution(wm, R_id, params, nlp_solver)

        if repaired
            for (a, connection) in wm.ref[:nw][n][:connection]
                network["pipes"][string(a)]["resistance"] = resistances[a][R_id[a]]
            end

            cvx = build_generic_model(network, CVXNLPWaterModel, WaterModels.post_cvx_hw)
            setsolver(cvx.model, nlp_solver)
            status = JuMP.solve(cvx.model)
            h = get_head_solution(cvx, n)

            # Set resistances appropriately.
            for (a, connection) in wm.ref[:nw][n][:connection]
                qp_sol = getvalue(cvx.var[:nw][n][:qp][a][1])
                qn_sol = getvalue(cvx.var[:nw][n][:qn][a][1])
                dir = (qp_sol - qn_sol) > 0.0 ? 1 : 0
                rneq = setdiff(1:length(resistances[a]), [R_id[a]])

                setsolutionvalue(cb, wm.var[:nw][n][:dir][a], dir)
                setsolutionvalue(cb, wm.var[:nw][n][:xr][a][R_id[a]], 1)

                for r in rneq
                    setsolutionvalue(cb, wm.var[:nw][n][:xr][a][r], 0)
                end

                dir_var = wm.var[:nw][n][:dir][a]
                L_a = connection["length"]

                # Add cuts.
                if qp_sol - qn_sol >= 0.0
                    qp = wm.var[:nw][n][:qp][a]
                    dhp = wm.var[:nw][n][:dhp][a]
                    lhs = compute_q_p_cut(dhp, qp, dir_var, qp_sol, resistances[a], R_id[a], L_a)
                    @usercut(cb, lhs <= 0.0)
                else
                    qn = wm.var[:nw][n][:qn][a]
                    dhn = wm.var[:nw][n][:dhn][a]
                    lhs = compute_q_n_cut(dhn, qn, dir_var, qn_sol, resistances[a], R_id[a], L_a)
                    @usercut(cb, lhs <= 0.0)
                end
            end

            addsolution(cb)
            params["obj_best"] = min(current_objective, params["obj_best"])
        end
    end

    return heuristic_cut_callback
end
