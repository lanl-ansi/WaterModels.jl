export lazy_cut_callback_generator

import MathProgBase

function lazy_cut_callback_generator(wm::GenericWaterModel,
                                     params::Dict{String, Any},
                                     nlp_solver::MathProgBase.AbstractMathProgSolver,
                                     n::Int = wm.cnw)
    network = deepcopy(wm.data)
    resistances = wm.ref[:nw][n][:resistance]
    connection_ids = collect(ids(wm, n, :connection))

    function lazy_cut_callback(cb::MathProgBase.MathProgCallbackData)
        # Set up variable arrays that will be used for cuts.
        xr_ones = Array{JuMP.Variable, 1}()
        xr_zeros = Array{JuMP.Variable, 1}()

        # Initialize the objective value.
        current_objective = 0.0
        R_id = Dict{Int, Int}(a => 0 for a in connection_ids)

        # Update resistances used throughout the network.
        for (a, connection) in wm.ref[:nw][n][:connection]
            xr_a = getvalue(wm.var[:nw][n][:xr][a])
            R_id[a] = findfirst(r -> isapprox(xr_a[r], 1.0, atol = 0.01), 1:length(xr_a))
            network["pipes"][string(a)]["resistance"] = resistances[a][R_id[a]]
            zero_indices = setdiff(1:length(xr_a), [R_id[a]])
            xr_ones = vcat(xr_ones, wm.var[:nw][n][:xr][a][R_id[a]])
            xr_zeros = vcat(xr_zeros, wm.var[:nw][n][:xr][a][zero_indices])
            L_a = wm.ref[:nw][n][:connection][a]["length"]
            current_objective += L_a * wm.ref[:nw][n][:resistance_cost][a][R_id[a]]
        end

        # Update objective values.
        params["obj_last"] = params["obj_curr"]
        params["obj_curr"] = current_objective

        # Solve the convex program.
        q, h = get_cvx_solution(wm, R_id, nlp_solver)
        qlb, qub, hlb, hub = check_solution_bounds(wm, q, h, R_id, n)
        println(all(values(qlb)), " ", all(values(qub)), " ", all(values(hlb)), " ", all(values(hub)))
        solution_is_feasible = all([all(collect(values(qlb))),
                                    all(collect(values(qub))),
                                    all(collect(values(hlb))),
                                    all(collect(values(hub)))])

        if !solution_is_feasible
            num_arcs = length(wm.ref[:nw][n][:connection])
            @lazyconstraint(cb, sum(xr_ones) - sum(xr_zeros) <= num_arcs - 1)
        else
            params["obj_best"] = min(current_objective, params["obj_best"])
        end
    end

    return lazy_cut_callback
end
