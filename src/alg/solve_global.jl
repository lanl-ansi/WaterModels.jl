export solve_global

const SolverType = MathProgBase.AbstractMathProgSolver

"""
Implements an algorithm to determine the globally optimal solution to a network
design problem. (This is a reproduction of Algorithm 1 in Raghunathan (2013).)
"""
function solve_global(network_path::String, problem_path::String,
                      nlp_solver::SolverType, mip_solver::SolverType)
    # Set the random seed.
    Random.seed!(1)

    # Load the required network data.
    network = WaterModels.parse_file(network_path)
    modifications = WaterModels.parse_file(problem_path)
    InfrastructureModels.update_data!(network, modifications)

    # Solve a relaxed version of the convex MINLP (Line 1).
    rnlp = build_generic_model(network, MINLPWaterModel, WaterModels.post_ne_hw)
    setsolver(rnlp.model, nlp_solver)
    rnlp_status = JuMP.solve(rnlp.model, relaxation = true)

    # Set initial points for the outer approximation (Lines 2 through 5).
    for (a, connection) in rnlp.ref[:nw][rnlp.cnw][:connection]
        # Get possible resistances for the arc.
        R_a = rnlp.ref[:nw][rnlp.cnw][:resistance][a]

        # Get outer approximation points for flow from i to j.
        q_p = getvalue(rnlp.var[:nw][rnlp.cnw][:qp][a])
        phi_max, r = findmax([R_a[r] * q_p[r]^(1.852) for r in 1:length(R_a)])

        # Get outer approximation points for flow from j to i.
        q_n = getvalue(rnlp.var[:nw][rnlp.cnw][:qn][a])
        phi_max, r = findmax([R_a[r] * q_n[r]^(1.852) for r in 1:length(R_a)])
    end

    # Make a dictionary for subalgorithm parameters (that are shared).
    obj_val = getobjectivevalue(rnlp.model)
    params = Dict{String, Any}("obj_last" => obj_val, "obj_curr" => obj_val,
                               "Beta_oa" => 5.0, "K_oa" => 1.0e-3,
                               "n" => Dict{Int, Int}(),
                               "max_repair_iters" => 50, "obj_best" => Inf,
                               "M_oa" => 5.0, "epsilon" => 1.0e-6)


    # Initialize the master MILP (mMILP) problem.
    mmilp = build_generic_model(network, MILPRWaterModel, WaterModels.post_ne_hw)

    ## Find an initial solution for the master MILP.
    #resistance_indices = find_initial_solution(mmilp, 20, 20, 0.85, nlp_solver)
    #objective_value = compute_objective(mmilp, resistance_indices)
    #println(resistance_indices, " ", objective_value)

    #return :Optimal

    #set_initial_solution(mmilp, best_resistance_indices, nlp_solver)

    # Set the solver for the problem and add the required callbacks.
    #user_cut_callback = user_cut_callback_generator(mmilp, params, mmilp.cnw)
    #addcutcallback(mmilp.model, user_cut_callback)
    lazy_cut_callback = lazy_cut_callback_generator(mmilp, params, nlp_solver)
    addlazycallback(mmilp.model, lazy_cut_callback)
    #heuristic_cut_callback = heuristic_cut_callback_generator(mmilp, params, nlp_solver, mmilp.cnw)
    #addheuristiccallback(mmilp.model, heuristic_cut_callback)

    # Solve the problem and return the status.
    setsolver(mmilp.model, mip_solver)
    return JuMP.solve(mmilp.model)
end
