export lazy_cut_callback_generator

import MathProgBase

function lazy_cut_callback_generator(wm::GenericWaterModel, params::Dict{String, Any},
                                     nlp_solver::MathProgBase.AbstractMathProgSolver,
                                     n::Int = wm.cnw)
    resistances = wm.ref[:nw][n][:resistance]
    network = deepcopy(wm.data)

    connection_ids = collect(ids(wm, n, :connection))
    R_id = Dict{Int, Int}(a => 1 for a in connection_ids)

    function lazy_cut_callback(cb::MathProgBase.MathProgCallbackData)
        # Set up variable arrays that will be used for cuts.
        xr_ones = Array{JuMP.Variable, 1}()
        xr_zeros = Array{JuMP.Variable, 1}()

        # Initialize the objective value.
        current_objective = 0.0

        # Update resistances used throughout the network.
        for (a, connection) in wm.ref[:nw][n][:connection]
            xr_a = getvalue(wm.var[:nw][n][:xr][a])
            R_id[a] = r = findfirst(r -> isapprox(xr_a[r], 1.0, atol = 0.01), 1:length(xr_a))
            network["pipes"][string(a)]["resistance"] = resistances[a][r]
            zero_indices = setdiff(1:length(xr_a), [r])
            xr_ones = vcat(xr_ones, wm.var[:nw][n][:xr][a][r])
            xr_zeros = vcat(xr_zeros, wm.var[:nw][n][:xr][a][zero_indices])
            L_a = wm.ref[:nw][n][:connection][a]["length"]
            current_objective += L_a * wm.ref[:nw][n][:resistance_cost][a][r]
        end

        # Update objective values.
        params["obj_last"] = params["obj_curr"]
        params["obj_curr"] = current_objective

        # Solve the convex program.
        sol_feasible, tmp = repair_solution(wm, R_id, params, nlp_solver,
                                            feasibility_check = true)

        if !sol_feasible
            num_arcs = length(wm.ref[:nw][n][:connection])
            @lazyconstraint(cb, sum(xr_ones) - sum(xr_zeros) <= num_arcs - 1)
        else
            params["obj_best"] = min(current_objective, params["obj_best"])
        end
    end

    return lazy_cut_callback
end
