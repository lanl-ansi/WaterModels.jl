using JuMP

export lazy_cut_callback_generator

import MathProgBase

function lazy_cut_callback_generator(wm::GenericWaterModel,
                                     params::Dict{String, Any},
                                     nlp_solver::MathProgBase.AbstractMathProgSolver,
                                     n::Int = wm.cnw)
    function lazy_cut_callback(cb::MathProgBase.MathProgCallbackData)
        # Set up variable arrays that will be used for cuts.
        xr_ones = Array{JuMP.Variable, 1}()
        xr_zeros = Array{JuMP.Variable, 1}()

        # Initialize the objective value.
        current_objective = 0.0
        resistances = wm.ref[:nw][n][:resistance]
        connection_ids = collect(ids(wm, n, :connection))
        resistance_indices = Dict{Int, Int}(a => 0 for a in connection_ids)

        # Update resistances used throughout the network.
        for (a, connection) in wm.ref[:nw][n][:connection]
            xr = getvalue(wm.var[:nw][n][:xr][a])
            r = findfirst(i -> isapprox(xr[i], 1.0; atol = 0.01), 1:length(xr))
            resistance_indices[a] = r

            zero_indices = setdiff(1:length(xr), [r])
            xr_ones = vcat(xr_ones, wm.var[:nw][n][:xr][a][r])
            xr_zeros = vcat(xr_zeros, wm.var[:nw][n][:xr][a][zero_indices])

            L_a = wm.ref[:nw][n][:connection][a]["length"]
            current_objective += L_a * wm.ref[:nw][n][:resistance_cost][a][r]
        end

        # Update objective values.
        params["obj_last"] = params["obj_curr"]
        params["obj_curr"] = current_objective

        # Solve the convex program.
        q, h = get_cvx_solution(wm, resistance_indices, nlp_solver)
        qlb, qub, hlb, hub = check_solution_bounds(wm, q, h, resistance_indices, n)
        solution_is_feasible = all([all(values(qlb)), all(values(qub)),
                                    all(values(hlb)), all(values(hub))])

        if !solution_is_feasible
            num_arcs = length(wm.ref[:nw][n][:connection])
            @lazyconstraint(cb, sum(xr_ones) - sum(xr_zeros) <= num_arcs - 1)
        else
            params["obj_best"] = min(current_objective, params["obj_best"])
        end
    end

    return lazy_cut_callback
end
