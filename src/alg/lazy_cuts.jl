export lazy_cut_callback_generator

import GLPK
import GLPKMathProgInterface
import MathProgBase

const SolverType = MathProgBase.AbstractMathProgSolver

function lazy_cut_callback_generator(wm::GenericWaterModel, params::Dict{String, Any}, nlp_solver::SolverType, n::Int = wm.cnw)
    resistances = wm.ref[:nw][n][:resistance]
    network = deepcopy(wm.data)

    function lazy_cut_callback(cb::MathProgBase.MathProgCallbackData)
        lp_solution = current_objective = nothing
        lp_solution = MathProgBase.cbgetmipsolution(cb)

        if typeof(cb) == GLPKMathProgInterface.GLPKInterfaceMIP.GLPKCallbackData
            lp_solution = MathProgBase.cbgetlpsolution(cb)
            current_problem = GLPK.ios_get_prob(cb.tree)
            current_objective = GLPK.get_obj_val(current_problem)
        end
        #elseif typeof(cb) == CPLEX.CplexLazyCallbackData
        #    lp_solution = MathProgBase.cbgetmipsolution(cb)
        #    current_objective = params["obj_curr"] #cbgetnodeobjval(cb)
        #end

        # Update objective values.
        params["obj_last"] = params["obj_curr"]
        params["obj_curr"] = 0.0 #current_objective

        # Set up variable arrays that will be used for cuts.
        xr_ones = Array{JuMP.Variable, 1}()
        xr_zeros = Array{JuMP.Variable, 1}()

        # Update resistances used throughout the network.
        for (a, connection) in wm.ref[:nw][n][:connection]
            xr_a = lp_solution[linearindex.(wm.var[:nw][n][:xr][a])]
            r = findfirst(r -> isapprox(xr_a[r], 1.0), 1:length(xr_a))
            network["pipes"][string(a)]["resistance"] = resistances[a][r]
            zero_indices = setdiff(1:length(xr_a), [r])
            xr_ones = vcat(xr_ones, wm.var[:nw][n][:xr][a][r])
            xr_zeros = vcat(xr_zeros, wm.var[:nw][n][:xr][a][zero_indices])
        end

        cvx = build_generic_model(network, MINLPWaterModel, WaterModels.post_cvx_hw)
        setsolver(cvx.model, nlp_solver)
        status = JuMP.solve(cvx.model, relaxation = true)

        if status != :LocalOptimal && status != :Optimal
            num_arcs = length(wm.ref[:nw][n][:connection])
            @lazyconstraint(cb, sum(xr_ones) - sum(xr_zeros) <= num_arcs - 1)
        elseif !WaterModels.solution_is_feasible(cvx, n)
            num_arcs = length(wm.ref[:nw][n][:connection])
            #@lazyconstraint(cb, sum(xr_ones) - sum(xr_zeros) <= num_arcs - 1)
        end
    end

    return lazy_cut_callback
end
