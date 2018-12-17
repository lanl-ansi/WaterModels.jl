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

    # Initialize the master MILP (mMILP) problem.
    mmilp = build_generic_model(network, MILPRWaterModel, WaterModels.post_ne_hw)

    # Find an initial solution for the master MILP.
    resistance_indices = find_initial_solution(mmilp, 50, 20, 0.85, nlp_solver)
    set_initial_solution(mmilp, resistance_indices, nlp_solver)
    best_objective = compute_objective(mmilp, resistance_indices)
    rnlp = build_generic_model(network, MINLPWaterModel, WaterModels.post_ne_hw)
    bound_tightening(mmilp, rnlp, best_objective, 2, nlp_solver)

    params = Dict{String, Any}("obj_last" => best_objective,
                               "obj_curr" => best_objective,
                               "obj_best" => best_objective,
                               "Beta_oa" => 5.0, "K_oa" => 1.0e-3,
                               "n" => Dict{Int, Int}(),
                               "max_repair_iters" => 50,
                               "M_oa" => 5.0, "epsilon" => 1.0e-6)

    println(resistance_indices, params["obj_best"])
    # Solve a relaxed version of the convex MINLP (Line 1).
    rnlp_cost = get_resistance_cost_expression(rnlp)
    @objective(rnlp.model, Min, rnlp_cost)
    setsolver(rnlp.model, nlp_solver)
    rnlp_status = JuMP.solve(rnlp.model, relaxation = true)

    # Set initial points for the outer approximation (Lines 2 through 5).
    for (a, connection) in rnlp.ref[:nw][rnlp.cnw][:connection]
        # Get possible resistances for the arc.
        R_a = mmilp.ref[:nw][mmilp.cnw][:resistance][a]
        L_a = mmilp.ref[:nw][mmilp.cnw][:connection][a]["length"]
        dir = mmilp.var[:nw][mmilp.cnw][:dir][a]

        # Add outer approximation cut for flow from i to j.
        q_p_sol = getvalue(rnlp.var[:nw][rnlp.cnw][:qp][a])
        phi_max_p, r_p = findmax([R_a[r] * q_p_sol[r]^(1.852) for r in 1:length(R_a)])
        qp = mmilp.var[:nw][mmilp.cnw][:qp][a]
        dhp = mmilp.var[:nw][mmilp.cnw][:dhp][a]
        lhs = compute_q_p_cut(dhp, qp, dir, q_p_sol[r_p], R_a, r_p, L_a)
        @constraint(mmilp.model, lhs <= 0.0)

        # Add outer approximation cut for flow from j to i.
        q_n_sol = getvalue(rnlp.var[:nw][rnlp.cnw][:qn][a])
        phi_max_n, r_n = findmax([R_a[r] * q_n_sol[r]^(1.852) for r in 1:length(R_a)])
        qn = mmilp.var[:nw][mmilp.cnw][:qn][a]
        dhn = mmilp.var[:nw][mmilp.cnw][:dhn][a]
        lhs = compute_q_n_cut(dhn, qn, dir, q_n_sol[r_n], R_a, r_n, L_a)
        @constraint(mmilp.model, lhs <= 0.0)
    end

    # Set the solver for the problem and add the required callbacks.
    user_cut_callback = user_cut_callback_generator(mmilp, params, nlp_solver)
    addcutcallback(mmilp.model, user_cut_callback)
    lazy_cut_callback = lazy_cut_callback_generator(mmilp, params, nlp_solver)
    addlazycallback(mmilp.model, lazy_cut_callback)
    heuristic_cut_callback = heuristic_cut_callback_generator(mmilp, params, nlp_solver, mmilp.cnw)
    addheuristiccallback(mmilp.model, heuristic_cut_callback)

    # Solve the problem and return the status.
    setsolver(mmilp.model, mip_solver)
    return JuMP.solve(mmilp.model)
end
