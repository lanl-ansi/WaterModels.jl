export lazy_cut_callback_generator

import GLPK

function lazy_cut_callback_generator(wm::GenericWaterModel, params::Dict{String, Any}, n::Int = wm.cnw)
    # Get indices associated with the MathProgBase model.
    arcs = collect(ids(wm, n, :connection))
    xr_indices = vcat([linearindex.(wm.var[:nw][n][:xr][a][:]) for a in arcs]...)
    network = deepcopy(wm.data)
    set_maximum_diameters(network)

    function lazy_cut_callback(cb::MathProgBase.MathProgCallbackData)
        lp_solution = MathProgBase.cbgetlpsolution(cb)
        xr_sol = lp_solution[xr_indices]

        current_problem = GLPK.ios_get_prob(cb.tree)
        current_objective = GLPK.get_obj_val(current_problem)
        println(current_objective)

        # Update objective values.
        params["obj_last"] = params["obj_curr"]
        params["obj_curr"] = current_objective

        # Update diameters in the network.
        #cvx = build_generic_model(network, MINLPWaterModel, WaterModels.post_wf_hw)

        #setsolver(cvx.model, IpoptSolver(linear_solver = "ma86"))
        #status = JuMP.solve(cvx.model, suppress_warnings = true)
        #WaterModels.solution_is_feasible(cvx)
    end

    return lazy_cut_callback
end
