export solve_global

const SolverType = MathProgBase.AbstractMathProgSolver

"""
Implements an algorithm to determine the globally optimal solution to a network
design problem. (This is a reproduction of Algorithm 1 in Raghunathan (2013).)
"""

function initialize_models(network_path::String, problem_path::String)
    # Load the required network data.
    network = WaterModels.parse_file(network_path)
    modifications = WaterModels.parse_file(problem_path)
    InfrastructureModels.update_data!(network, modifications)

    # Initialize the master MILP (mMILP) and relaxed NLP (RNLP) problems.
    mmilp = build_generic_model(network, MILPRWaterModel, WaterModels.post_ne_hw)
    rnlp = build_generic_model(network, MINLPWaterModel, WaterModels.post_ne_hw)

    # Return the models.
    return mmilp, rnlp
end

function solve_global(network_path::String, problem_path::String,
                      nlp_solver::SolverType, mip_solver::SolverType)
    # Set the random seed.
    Random.seed!(1)

    # Initialize the master MILP and relaxed NLP models.
    (mmilp, rnlp), initialize_models_time = @timed initialize_models(network_path, problem_path)

    # Eliminate unnecessary binary variables.
    nothing, eliminate_variables_time = @timed eliminate_variables(mmilp, rnlp, nlp_solver)

    # Find an initial solution for the master MILP.
    resistance_indices, find_initial_solution_time = @timed find_initial_solution(mmilp, 250, 100, 0.90, nlp_solver)

    # Set the initial solution.
    nothing, set_initial_solution_time = @timed set_initial_solution(mmilp, resistance_indices, nlp_solver)
    best_objective = compute_objective(mmilp, resistance_indices)
    objective_initial = best_objective

    # Initialize the outer approximation.
    nothing, add_outer_approximation_time = @timed add_outer_approximation(mmilp, rnlp, resistance_indices, nlp_solver)

    # Initialize solving parameters.
    params = Dict{String, Any}("obj_last" => best_objective,
                               "obj_curr" => best_objective,
                               "obj_best" => best_objective,
                               "Beta_oa" => 5.0, "K_oa" => 1.0e-3,
                               "max_repair_iters" => 50,
                               "M_oa" => 5.0, "epsilon" => 1.0e-6)

    # Set the solver for the problem and add the required callbacks.
    user_cut_callback = user_cut_callback_generator(mmilp, params, nlp_solver)
    addcutcallback(mmilp.model, user_cut_callback)
    lazy_cut_callback = lazy_cut_callback_generator(mmilp, params, nlp_solver)
    addlazycallback(mmilp.model, lazy_cut_callback)
    heuristic_cut_callback = heuristic_cut_callback_generator(mmilp, params, nlp_solver)
    addheuristiccallback(mmilp.model, heuristic_cut_callback)

    # Solve the problem and return the status.
    setsolver(mmilp.model, mip_solver)
    solve_status, solve_time = @timed JuMP.solve(mmilp.model)

    # Get objective value information.
    objective_value = getobjectivevalue(mmilp.model)
    objective_bound = getobjbound(mmilp.model)
    objective_gap = abs(objective_bound - objective_value) / abs(objective_value)

    # Get the total time.
    total_time = initialize_models_time + find_initial_solution_time +
                 set_initial_solution_time + add_outer_approximation_time +
                 solve_time

    println("--------------------------- STATISTICS ---------------------------")
    println("initialize_models_time       = ", initialize_models_time)
    println("eliminate_variables_time     = ", eliminate_variables_time)
    println("find_initial_solution_time   = ", find_initial_solution_time)
    println("set_initial_solution_time    = ", set_initial_solution_time)
    println("add_outer_approximation_time = ", add_outer_approximation_time)
    println("solve_time                   = ", solve_time)
    println("total_time                   = ", total_time)
    println("solve_status                 = ", solve_status)
    println("objective_initial            = ", objective_initial)
    println("objective_value              = ", objective_value)
    println("objective_bound              = ", objective_bound)
    println("objective_gap                = ", objective_gap)
    println("------------------------------------------------------------------")

    return solve_status
end
