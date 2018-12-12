export repair_solution

import MathProgBase

"""
Implements an algorithm to repair a set of resistance choices whose network
analysis solution violates boundss on flow or potentials (Algorithm 3 in
Raghunathan (2013)). Returns a repaired set of resistance choices.
"""
function repair_solution(network::Dict{String, Any}, R::Dict{Int, Array{Float64, 1}}, R_id::Dict{Int, Int}, params::Dict{String, Any}, nlp_solver::MathProgBase.AbstractMathProgSolver)
    num_iterations = 1
    progress = true
    repaired = false

    while num_iterations <= params["max_repair_iters"] && !repaired && progress
        ## Update resistances.
        #for a in collect(keys(R))
        #    network["pipes"][string(a)]["resistance"] = R[a][R_id[a]]
        #end

        #cvx = build_generic_model(network, CVXNLPWaterModel, WaterModels.post_cvx_hw)
        #setsolver(cvx.model, nlp_solver)
        #status = JuMP.solve(cvx.model, relaxation = true, suppress_warnings = true)

        #if
        #end

        num_iterations += 1
    end

    #println(network)
    #println(resistances)
    #println(params)
    #r = resistances

    #while num_iterations <= max_iterations && !repaired
    #    num_iterations += 1
    #end
end
