export lazy_cut_callback_generator

import GLPK

function lazy_cut_callback_generator(wm::GenericWaterModel, params::Dict{String, Any}, n::Int = wm.cnw)
    # Get indices associated with the MathProgBase model.
    arcs = collect(ids(wm, n, :connection))
    dir_indices = linearindex.(wm.var[:nw][n][:dir][:])
    xr_indices = vcat([linearindex.(wm.var[:nw][n][:xr][a][:]) for a in arcs]...)

    function lazy_cut_callback(cb::MathProgBase.MathProgCallbackData)
        lp_solution = MathProgBase.cbgetlpsolution(cb)
        dir_sol = lp_solution[dir_indices]
        xr_sol = lp_solution[xr_indices]

        # TODO: The method for finding the node level should be solver-agnostic.
        current_node = GLPK.ios_curr_node(cb.tree)
        previous_node = GLPK.ios_prev_node(cb.tree, current_node)
        node_depth = GLPK.ios_node_level(cb.tree, current_node)
        current_problem = GLPK.ios_get_prob(cb.tree)
        current_objective = GLPK.get_obj_val(current_problem)

        # Update objective values.
        params["obj_last"] = params["obj_curr"]
        params["obj_curr"] = current_objective
    end

    return lazy_cut_callback
end
