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
    mmilp = build_generic_model(network, MILPRWaterModel, WaterModels.post_ne_hw_segmented)
    rnlp = build_generic_model(network, MINLPWaterModel, WaterModels.post_ne_hw)

    ## Eliminate unnecessary binary variables.
    #eliminate_variables(mmilp, rnlp, nlp_solver)

    # Find an initial solution for the master MILP.
    resistance_indices = find_initial_solution(mmilp, 100, 50, 0.80, nlp_solver)
    set_initial_solution(mmilp, 3, resistance_indices, nlp_solver)
    best_objective = compute_objective(mmilp, resistance_indices)
    println("Starting objective value: ", best_objective)

    get_partitioning(mmilp, resistance_indices, nlp_solver)
    return :Optimal

    # Tighten bounds.
    #bound_tightening(mmilp, rnlp, best_objective, 2, nlp_solver)

    #params = Dict{String, Any}("obj_last" => best_objective,
    #                           "obj_curr" => best_objective,
    #                           "obj_best" => best_objective,
    #                           "Beta_oa" => 5.0, "K_oa" => 1.0e-3,
    #                           "max_repair_iters" => 50,
    #                           "M_oa" => 5.0, "epsilon" => 1.0e-6)

    ## Solve a relaxed version of the convex MINLP (Line 1).
    #rnlp_cost = get_resistance_cost_expression(rnlp)
    #@objective(rnlp.model, Min, rnlp_cost)
    #setsolver(rnlp.model, nlp_solver)
    #rnlp_status = JuMP.solve(rnlp.model, relaxation = true)
    #relaxed_objective_value = getobjectivevalue(rnlp.model)
    #println("Relaxed NLP objective value: ", relaxed_objective_value)

    ## Set initial points for the outer approximation (Lines 2 through 5).
    #for (a, connection) in mmilp.ref[:nw][mmilp.cnw][:connection]
    #    # Get possible resistances for the arc.
    #    R_a = mmilp.ref[:nw][mmilp.cnw][:resistance][a]
    #    L_a = mmilp.ref[:nw][mmilp.cnw][:connection][a]["length"]
    #    dir = mmilp.var[:nw][mmilp.cnw][:dir][a]

    #    # Add outer approximation cut for flow from i to j.
    #    q_p_sol = max.(getvalue(rnlp.var[:nw][rnlp.cnw][:qp][a]), 0.0)
    #    phi_max_p, r_p = findmax([R_a[r] * q_p_sol[r]^(1.852) for r in 1:length(R_a)])
    #    qp = mmilp.var[:nw][mmilp.cnw][:qp][a]
    #    dhp = mmilp.var[:nw][mmilp.cnw][:dhp][a]
    #    lhs = compute_q_p_cut(dhp, qp, dir, q_p_sol[r_p], R_a, r_p, L_a)
    #    @constraint(mmilp.model, lhs <= 0.0)

    #    # Add outer approximation cut for flow from j to i.
    #    q_n_sol = max.(getvalue(rnlp.var[:nw][rnlp.cnw][:qn][a]), 0.0)
    #    phi_max_n, r_n = findmax([R_a[r] * q_n_sol[r]^(1.852) for r in 1:length(R_a)])
    #    qn = mmilp.var[:nw][mmilp.cnw][:qn][a]
    #    dhn = mmilp.var[:nw][mmilp.cnw][:dhn][a]
    #    lhs = compute_q_n_cut(dhn, qn, dir, q_n_sol[r_n], R_a, r_n, L_a)
    #    @constraint(mmilp.model, lhs <= 0.0)
    #end

    ## Set the solver for the problem and add the required callbacks.
    #user_cut_callback = user_cut_callback_generator(mmilp, params, nlp_solver)
    #addcutcallback(mmilp.model, user_cut_callback)
    #lazy_cut_callback = lazy_cut_callback_generator(mmilp, params, nlp_solver)
    #addlazycallback(mmilp.model, lazy_cut_callback)
    #heuristic_cut_callback = heuristic_cut_callback_generator(mmilp, params, nlp_solver, 3, mmilp.cnw)
    #addheuristiccallback(mmilp.model, heuristic_cut_callback)

    ## Solve the problem and return the status.
    #setsolver(mmilp.model, mip_solver)
    #status = JuMP.solve(mmilp.model)
    #final_objective_value = getobjectivevalue(mmilp.model)
    #println("Final objective value: ", final_objective_value)
    #return status
end
