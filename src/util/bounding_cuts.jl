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

function generate_breakpoints(x::JuMP.Variable, x_c::Float64, num_breakpoints::Int)
    num_breakpoints_left = floor(num_breakpoints / 2)
    num_breakpoints_right = floor(num_breakpoints / 2)

    if isodd(num_breakpoints)
        num_breakpoints_left += 1
    end

    x_lb = getlowerbound(x)
    x_ub = getupperbound(x)

    if x_ub - x_lb <= 0.0
        return zeros(Float64, num_breakpoints)
    end

    x_c = x_c > x_lb + 1.0e-6 ? x_c : 0.5 * (x_ub - x_lb)
    x_c = x_c < x_ub - 1.0e-6 ? x_c : 0.5 * (x_ub - x_lb)

    left_width = x_c - x_lb
    initial_width = 0.001 < left_width ? 0.001 : left_width / 100.0
    exponent = log(left_width / initial_width) / (num_breakpoints_left - 1)
    c = initial_width * exp(-log(left_width / initial_width) / (num_breakpoints_left - 1))
    x_l = x_c .- [c * exp(i*exponent) for i in 1:num_breakpoints_left]

    right_width = x_ub - x_c
    initial_width = 0.001 < right_width ? 0.001 : right_width / 100.0
    exponent = log(right_width / initial_width) / (num_breakpoints_right - 1)
    c = initial_width * exp(-log(right_width / initial_width) / (num_breakpoints_right - 1))
    x_r = x_c .+ [c * exp(i*exponent) for i in 1:num_breakpoints_right]

    breakpoints = sort(vcat(x_l, x_r))
    zero_indices = findall(k -> breakpoints[k] < 1.0e-6, 1:length(breakpoints))
    setindex!(breakpoints, zeros(size(zero_indices, 1)), zero_indices)
    return breakpoints
end

function add_upper_approximation(wm::GenericWaterModel,
                                 rwm::GenericWaterModel,
                                 resistance_indices::Dict{Int, Int},
                                 solver::MathProgBase.AbstractMathProgSolver,
                                 n::Int = wm.cnw)
    q, h = get_cvx_solution(wm, resistance_indices, solver)
    num_breakpoints = 4

    # Set initial points for the outer approximation (Lines 2 through 5).
    for (a, connection) in wm.ref[:nw][n][:connection]
        # Get possible resistances for the arc.
        R_a = wm.ref[:nw][n][:resistance][a]
        dir = wm.var[:nw][n][:dir][a]
        r_hat = resistance_indices[a]

        # Create partition variables and constrain the number of partitions.
        x_par = @variable(wm.model, [k in 1:num_breakpoints-1], category = :Bin, start = 0)

        con_x_par = @constraint(wm.model, sum(x_par) == 1)

        lambda_p = @variable(wm.model, [r in 1:length(R_a), k in 1:num_breakpoints],
                             category = :Cont, lowerbound = 0.0,
                             upperbound = 1.0, start = 0.0)

        lambda_n = @variable(wm.model, [r in 1:length(R_a), k in 1:num_breakpoints],
                             category = :Cont, lowerbound = 0.0,
                             upperbound = 1.0, start = 0.0)

        f_hat_p = zero(AffExpr)
        f_hat_n = zero(AffExpr)

        for r in 1:length(R_a)
            x_res = wm.var[:nw][n][:xr][a][r]

            con_lambda_1 = @constraint(wm.model, sum(lambda_p[r, :]) + sum(lambda_n[r, :]) == x_res)
            con_lambda_n_1 = @constraint(wm.model, sum(lambda_n[r, :]) <= (1 - dir))
            con_lambda_p_1 = @constraint(wm.model, sum(lambda_p[r, :]) <= dir)
            con_lambda_n_2 = @constraint(wm.model, lambda_n[r, 1] <= x_par[1])
            con_lambda_p_2 = @constraint(wm.model, lambda_p[r, 1] <= x_par[1])

            for k in 2:num_breakpoints-1
                @constraint(wm.model, lambda_n[r, k] <= x_par[k-1] + x_par[k])
                @constraint(wm.model, lambda_p[r, k] <= x_par[k-1] + x_par[k])
            end

            con_lambda_n_3 = @constraint(wm.model, lambda_n[r, end] <= x_par[end])
            con_lambda_p_3 = @constraint(wm.model, lambda_p[r, end] <= x_par[end])

            qp_c = max(0.0, getvalue(rwm.var[:nw][n][:qp][a][r]))
            qp = generate_breakpoints(wm.var[:nw][n][:qp][a][r], qp_c, num_breakpoints)
            fp = [R_a[r] * qp[k]^1.852 for k in 1:num_breakpoints]
            f_hat_p += AffExpr(lambda_p[r, :], fp, 0.0)

            qn_c = max(0.0, getvalue(rwm.var[:nw][n][:qn][a][r]))
            qn = generate_breakpoints(wm.var[:nw][n][:qn][a][r], qn_c, num_breakpoints)
            fn = [R_a[r] * qn[k]^1.852 for k in 1:num_breakpoints]
            f_hat_n += AffExpr(lambda_n[r, :], fn, 0.0)

            for k in 1:num_breakpoints-1
                if q[a] > 0.0 && r == r_hat
                    if q[a] >= qp[k] && q[a] <= qp[k+1]
                        setvalue(x_par[k], 1)
                        lambda_val = (q[a] - qp[k]) / (qp[k+1] - qp[k])
                        setvalue(lambda_p[r, k+1], lambda_val)
                        setvalue(lambda_p[r, k], 1.0 - lambda_val)
                    end
                elseif q[a] <= 0.0 && r == r_hat
                    if -q[a] >= qn[k] && -q[a] <= qn[k+1]
                        setvalue(x_par[k], 1)
                        lambda_val = (-q[a] - qn[k]) / (qn[k+1] - qn[k])
                        setvalue(lambda_n[r, k+1], lambda_val)
                        setvalue(lambda_n[r, k], 1.0 - lambda_val)
                    end
                end
            end
        end

        dhp = wm.var[:nw][n][:dhp][a]
        dhn = wm.var[:nw][n][:dhn][a]
        L_a = wm.ref[:nw][n][:connection][a]["length"]

        @constraint(wm.model, dhp / L_a <= f_hat_p)
        @constraint(wm.model, dhn / L_a <= f_hat_n)
    end
end
