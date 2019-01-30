function add_outer_approximation(wm::GenericWaterModel,
                                 rwm::GenericWaterModel,
                                 resistance_indices::Dict{Int, Int},
                                 solver::MathProgBase.AbstractMathProgSolver,
                                 n::Int = wm.cnw)
    # Solve a relaxed version of the convex MINLP (Line 1).
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

function add_upper_approximation(wm::GenericWaterModel,
                                 resistance_indices::Dict{Int, Int},
                                 solver::MathProgBase.AbstractMathProgSolver,
                                 n::Int = wm.cnw)
    q, h = get_cvx_solution(wm, resistance_indices, solver)

    # Set initial points for the outer approximation (Lines 2 through 5).
    for (a, connection) in wm.ref[:nw][n][:connection]
        # Get possible resistances for the arc.
        R_a = wm.ref[:nw][n][:resistance][a]
        L_a = wm.ref[:nw][n][:connection][a]["length"]
        dir = wm.var[:nw][n][:dir][a]
        r_hat = resistance_indices[a]
        x_res = wm.var[:nw][n][:xr][a][r_hat]

        # Select a set of resistances where upper approximations will be applied.
        r_to_bound = intersect(1:length(R_a), [r_hat])
        r_ignored = []

        if q[a] > 0.0
            bound_expr = zero(AffExpr)

            for r in r_to_bound
                q_sol = compute_q_tilde(q[a], R_a[r_hat], R_a[r])

                # Set up partitioning points.
                q_1 = getlowerbound(wm.var[:nw][n][:qp][a][r])
                q_4 = getupperbound(wm.var[:nw][n][:qp][a][r])

                if r == r_hat
                    q_2 = max(q_1, q_sol - (q_4 - q_1) / 1000.0)
                    q_3 = min(q_4, q_sol + (q_4 - q_1) / 1000.0)
                else
                    q_2 = max(q_1, q_sol - (q_4 - q_1) / 100.0)
                    q_3 = min(q_4, q_sol + (q_4 - q_1) / 100.0)
                end

                if (q_1 == q_2 || q_3 == q_4)
                    r_ignored = vcat(r_ignored, r)
                    continue
                end

                q_par = @variable(wm.model, [k in 1:3], category = :Cont, lowerbound = 0.0, start = 0.0)
                x_par = @variable(wm.model, [k in 1:3], category = :Bin, start = 0)

                if r == r_hat
                    setvalue(x_par[2], 1)
                    setvalue(q_par[2], q[a])
                end

                qp_lbs = [q_1, q_2, q_3]
                qp_ubs = [q_2, q_3, q_4]
                setupperbound(q_par[1], q_2)
                setupperbound(q_par[2], q_3)
                setupperbound(q_par[3], q_4)

                con_1 = @constraint(wm.model, sum(q_par) == wm.var[:nw][n][:qp][a][r])
                con_2 = @constraint(wm.model, q_par[1] >= x_par[1] * q_1)
                con_3 = @constraint(wm.model, q_par[1] <= x_par[1] * q_2)
                con_4 = @constraint(wm.model, q_par[2] >= x_par[2] * q_2)
                con_5 = @constraint(wm.model, q_par[2] <= x_par[2] * q_3)
                con_6 = @constraint(wm.model, q_par[3] >= x_par[3] * q_3)
                con_7 = @constraint(wm.model, q_par[3] <= x_par[3] * q_4)

                dfp = R_a[r] .* [qp_ubs[k]^1.852 - qp_lbs[k]^1.852 for k in 1:3]
                dqp = [qp_ubs[k] - qp_lbs[k] for k in 1:3]
                qp_slopes = [dfp[k] / dqp[k] for k in 1:3]
                nan_indices = findall(isnan, qp_slopes)
                setindex!(qp_slopes, zeros(size(nan_indices)), nan_indices)
                qp_constants = [R_a[r] * qp_lbs[k]^1.852 - qp_lbs[k] * qp_slopes[k] for k in 1:3]

                term_1 = AffExpr(q_par, qp_slopes, 0.0)
                term_2 = AffExpr(x_par, qp_constants, 0.0)
                bound_expr += term_1 + term_2
            end

            dhp = wm.var[:nw][n][:dhp][a]
            r_rest = vcat(setdiff(1:length(R_a), r_to_bound), r_ignored)
            qp_ubs = [getupperbound(wm.var[:nw][n][:qp][a][r]) for r in 1:length(R_a)]
            slopes_p = [R_a[r]*qp_ubs[r]^(0.852) for r in 1:length(R_a)]
            remainder = AffExpr(wm.var[:nw][n][:qp][a][r_rest], slopes_p[r_rest], 0.0)
            con_8 = @constraint(wm.model, dhp / L_a <= bound_expr + remainder)
        elseif q[a] < 0.0
            bound_expr = zero(AffExpr)

            for r in r_to_bound
                q_sol = compute_q_tilde(-q[a], R_a[r_hat], R_a[r])

                # Set up partitioning points.
                q_1 = getlowerbound(wm.var[:nw][n][:qn][a][r])
                q_4 = getupperbound(wm.var[:nw][n][:qn][a][r])

                if r == r_hat
                    q_2 = max(q_1, q_sol - (q_4 - q_1) / 1000.0)
                    q_3 = min(q_4, q_sol + (q_4 - q_1) / 1000.0)
                else
                    q_2 = max(q_1, q_sol - (q_4 - q_1) / 100.0)
                    q_3 = min(q_4, q_sol + (q_4 - q_1) / 100.0)
                end

                if (q_1 == q_2 || q_3 == q_4)
                    r_ignored = vcat(r_ignored, r)
                    continue
                end

                q_par = @variable(wm.model, [k in 1:3], category = :Cont, lowerbound = 0.0, start = 0.0)
                x_par = @variable(wm.model, [k in 1:3], category = :Bin, start = 0)

                if r == r_hat
                    setvalue(x_par[2], 1)
                    setvalue(q_par[2], -q[a])
                end

                qn_lbs = [q_1, q_2, q_3]
                qn_ubs = [q_2, q_3, q_4]
                setupperbound(q_par[1], q_2)
                setupperbound(q_par[2], q_3)
                setupperbound(q_par[3], q_4)

                con_1 = @constraint(wm.model, sum(q_par) == wm.var[:nw][n][:qn][a][r])
                con_2 = @constraint(wm.model, q_par[1] >= x_par[1] * q_1)
                con_3 = @constraint(wm.model, q_par[1] <= x_par[1] * q_2)
                con_4 = @constraint(wm.model, q_par[2] >= x_par[2] * q_2)
                con_5 = @constraint(wm.model, q_par[2] <= x_par[2] * q_3)
                con_6 = @constraint(wm.model, q_par[3] >= x_par[3] * q_3)
                con_7 = @constraint(wm.model, q_par[3] <= x_par[3] * q_4)

                dfn = R_a[r] .* [qn_ubs[k]^1.852 - qn_lbs[k]^1.852 for k in 1:3]
                dqn = [qn_ubs[k] - qn_lbs[k] for k in 1:3]
                qn_slopes = [dfn[k] / dqn[k] for k in 1:3]
                nan_indices = findall(isnan, qn_slopes)
                setindex!(qn_slopes, zeros(size(nan_indices)), nan_indices)
                qn_constants = [R_a[r] * qn_lbs[k]^1.852 - qn_lbs[k] * qn_slopes[k] for k in 1:3]

                term_1 = AffExpr(q_par, qn_slopes, 0.0)
                term_2 = AffExpr(x_par, qn_constants, 0.0)
                bound_expr += term_1 + term_2
            end

            dhn = wm.var[:nw][n][:dhn][a]
            r_rest = vcat(setdiff(1:length(R_a), r_to_bound), r_ignored)
            qn_ubs = [getupperbound(wm.var[:nw][n][:qn][a][r]) for r in 1:length(R_a)]
            slopes_n = [R_a[r]*qn_ubs[r]^(0.852) for r in 1:length(R_a)]
            remainder = AffExpr(wm.var[:nw][n][:qn][a][r_rest], slopes_n[r_rest], 0.0)
            con_8 = @constraint(wm.model, dhn / L_a <= bound_expr + remainder)
        end
    end
end
