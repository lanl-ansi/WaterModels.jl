function add_outer_approximation(wm::GenericWaterModel,
                                 rwm::GenericWaterModel,
                                 resistance_indices::Dict{Int, Int},
                                 solver::MathProgBase.AbstractMathProgSolver,
                                 n::Int = wm.cnw)
    # Solve a relaxed version of the MICP (Line 1).
    rwm_cost = get_resistance_cost_expression(rwm)
    @objective(rwm.model, Min, rwm_cost)
    setsolver(rwm.model, solver)
    rwm_status = JuMP.solve(rwm.model, relaxation = true)

    # Set initial points for the outer approximation (Lines 2 through 5).
    for (a, connection) in wm.ref[:nw][n][:connection]
        # Get possible resistances for the arc.
        R_a = wm.ref[:nw][n][:resistance][a]
        L_a = wm.ref[:nw][n][:connection][a]["length"]
        dir = wm.var[:nw][n][:dir][a]

        # Add outer approximation cut for flow from i to j.
        q_p_sol = getvalue(rwm.var[:nw][n][:qp][a])
        phi_max_p, r_p = findmax([R_a[r] * q_p_sol[r]^(1.852) for r in 1:length(R_a)])
        qp = wm.var[:nw][n][:qp][a]
        dhp = wm.var[:nw][n][:dhp][a]
        lhs = compute_q_p_cut(dhp, qp, dir, q_p_sol[r_p], R_a, r_p, L_a)
        @constraint(wm.model, lhs <= 0.0)

        # Add outer approximation cut for flow from j to i.
        q_n_sol = getvalue(rwm.var[:nw][n][:qn][a])
        phi_max_n, r_n = findmax([R_a[r] * q_n_sol[r]^(1.852) for r in 1:length(R_a)])
        qn = wm.var[:nw][n][:qn][a]
        dhn = wm.var[:nw][n][:dhn][a]
        lhs = compute_q_n_cut(dhn, qn, dir, q_n_sol[r_n], R_a, r_n, L_a)
        @constraint(wm.model, lhs <= 0.0)
    end
end
