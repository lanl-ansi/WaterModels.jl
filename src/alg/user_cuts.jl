export user_cut_callback_generator

import GLPK
import Random

function user_cut_callback_generator(wm::GenericWaterModel, params::Dict{String, Any}, n::Int = wm.cnw)
    # Get indices associated with the MathProgBase model.
    arcs = collect(ids(wm, n, :connection))
    dir_indices = linearindex.(wm.var[:nw][n][:dir][:])
    xr_indices = vcat([linearindex.(wm.var[:nw][n][:xr][a][:]) for a in arcs]...)

    function user_cut_callback(cb::MathProgBase.MathProgCallbackData)
        lp_solution = MathProgBase.cbgetlpsolution(cb)
        dir_sol = lp_solution[dir_indices]
        xr_sol = lp_solution[xr_indices]

        # TODO: The method for finding these data should be solver-agnostic.
        current_node = GLPK.ios_curr_node(cb.tree)
        previous_node = GLPK.ios_prev_node(cb.tree, current_node)
        current_problem = GLPK.ios_get_prob(cb.tree)
        current_objective = GLPK.get_obj_val(current_problem)
        d = convert(Float64, GLPK.ios_node_level(cb.tree, current_node))

        # Update objective values.
        params["obj_last"] = params["obj_curr"]
        params["obj_curr"] = current_objective

        # Conditions for adding outer approximations.
        depth_satisfied = Random.rand() <= params["Beta_oa"] * 2^(-d)
        obj_change_rel = (params["obj_curr"] - params["obj_last"]) / params["obj_last"]
        obj_improved = obj_change_rel >= params["K_oa"]
       
        # Initialize the number of cut rounds added per node.
        if !haskey(params["n"], current_node)
            params["n"][current_node] = 0
        end

        # Check satisfaction of the number of rounds.
        num_rounds_satisfied = params["n"][current_node] <= params["M_oa"]

        if depth_satisfied && obj_improved && num_rounds_satisfied
            for (relative_index, a) in enumerate(arcs)
                if dir_sol[relative_index] >= 0.5
                    #if violation > params["epsilon"]
                    #    # params["oa_p"]
                    #end
                else
                    #if violation > params["epsilon"]
                    #    # params["oa_n"]
                    #end
                end
            end

            params["n"][current_node] += 1
        end
    end

    return user_cut_callback
end
