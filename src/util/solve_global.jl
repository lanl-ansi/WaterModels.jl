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

    # Initialize the master MILP (mMILP) and relaxed NLP (RNLP) problems.
    mmilp = build_generic_model(network, MILPRWaterModel, WaterModels.post_ne_hw)
    rnlp = build_generic_model(network, MINLPWaterModel, WaterModels.post_ne_hw)

    ## Eliminate unnecessary binary variables.
    #eliminate_variables(mmilp, rnlp, nlp_solver)

    # Find an initial solution for the master MILP.
    resistance_indices = find_initial_solution(mmilp, 100, 50, 0.80, nlp_solver)
    set_initial_solution(mmilp, resistance_indices, nlp_solver)
    best_objective = compute_objective(mmilp, resistance_indices)
    println("Starting objective value: ", best_objective)

    ## Tighten bounds.
    #bound_tightening(mmilp, rnlp, best_objective, 2, nlp_solver)

    params = Dict{String, Any}("obj_last" => best_objective,
                               "obj_curr" => best_objective,
                               "obj_best" => best_objective,
                               "Beta_oa" => 5.0, "K_oa" => 1.0e-3,
                               "max_repair_iters" => 50,
                               "M_oa" => 5.0, "epsilon" => 1.0e-6)

    # Initialize the outer approximation.
    add_outer_approximation(mmilp, rnlp, resistance_indices, nlp_solver)
    add_upper_approximation(mmilp, resistance_indices, nlp_solver)

    # Set the solver for the problem and add the required callbacks.
    user_cut_callback = user_cut_callback_generator(mmilp, params, nlp_solver)
    addcutcallback(mmilp.model, user_cut_callback)
    lazy_cut_callback = lazy_cut_callback_generator(mmilp, params, nlp_solver)
    addlazycallback(mmilp.model, lazy_cut_callback)
    ##heuristic_cut_callback = heuristic_cut_callback_generator(mmilp, params, nlp_solver, 3, mmilp.cnw)
    ##addheuristiccallback(mmilp.model, heuristic_cut_callback)

    # Solve the problem and return the status.
    setsolver(mmilp.model, mip_solver)
    status = JuMP.solve(mmilp.model)
    final_objective_value = getobjectivevalue(mmilp.model)
    println("Final objective value: ", final_objective_value)
    return status
end
