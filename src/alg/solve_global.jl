export solve_global

const SolverType = MathProgBase.AbstractMathProgSolver

"""
Implements an algorithm to determine the globally optimal solution to a network
design problem. (This is a reproduction of Algorithm 1 in Raghunathan (2013).)
"""
function solve_global(network_path::String, problem_path::String, nlp_solver::SolverType, mip_solver::SolverType)
    # Load the required network data.
    network = WaterModels.parse_file(network_path)
    modifications = WaterModels.parse_file(problem_path)
    InfrastructureModels.update_data!(network, modifications)

    # Solve a relaxed version of the convex MINLP (Line 1).
    rnlp = build_generic_model(network, MINLPWaterModel, WaterModels.post_ne_hw)
    setsolver(rnlp.model, nlp_solver)
    rnlp_status = JuMP.solve(rnlp.model, relaxation = true)

    # Get the dictionary of resistances for the problem.
    R = calc_resistances_hw(rnlp, rnlp.cnw)

    # Initialize the set of points at which outer approximations will be made.
    oa_p = [Set{Tuple{Float64, Int}}() for a in collect(ids(rnlp, :connection))]
    oa_n = [Set{Tuple{Float64, Int}}() for a in collect(ids(rnlp, :connection))]

    # Set initial points for the outer approximation (Lines 2 through 5).
    for (a, connection) in rnlp.ref[:nw][rnlp.cnw][:connection]
        # Get outer approximation points for flow from i to j.
        q_p = getvalue(rnlp.var[:nw][rnlp.cnw][:qp][a])
        phi_max, r = findmax([R[a][r] * q_p[r]^(1.852) for r in 1:length(q_p)])
        push!(oa_p[a], (q_p[r], r))

        # Get outer approximation points for flow from j to i.
        q_n = getvalue(rnlp.var[:nw][rnlp.cnw][:qn][a])
        phi_max, r = findmax([R[a][r] * q_n[r]^(1.852) for r in 1:length(q_n)])
        push!(oa_n[a], (q_n[r], r))
    end

    # Make a dictionary for subalgorithm parameters (that are shared).
    obj_val = getobjectivevalue(rnlp.model)
    params = Dict{String, Any}("obj_last" => obj_val, "obj_curr" => obj_val,
                               "oa_p" => oa_p, "oa_n" => oa_n, "Beta_oa" => 5.0,
                               "K_oa" => 1.0e-3, "n" => Dict{Int, Int}(),
                               "max_repair_iters" => 50, "M_oa" => 5.0,
                               "epsilon" => 1.0e-6)

    # Initialize the master MILP (mMILP) problem.
    mmilp = build_generic_model(network, MILPRWaterModel, WaterModels.post_ne_hw)

    # Set the solver for the problem and add the required callbacks.
    setsolver(mmilp.model, mip_solver)
    user_cut_callback = user_cut_callback_generator(mmilp, params, mmilp.cnw)
    addcutcallback(mmilp.model, user_cut_callback)
    lazy_cut_callback = lazy_cut_callback_generator(mmilp, params, nlp_solver, mmilp.cnw)
    addlazycallback(mmilp.model, lazy_cut_callback)

    # Solve the problem and return the status.
    return JuMP.solve(mmilp.model, relaxation = false)
end
