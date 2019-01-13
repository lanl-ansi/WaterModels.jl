export user_cut_callback_generator, compute_q_tilde, compute_q_p_cut, compute_q_n_cut

import Random

function compute_q_tilde(q_hat::Float64, r_hat::Float64, r::Float64)
    b = r_hat * (head_loss_hw_prime(q_hat) * q_hat - head_loss_hw_func(q_hat))
    return 1.0903341484 * (b / r)^0.53995680345
end

function compute_q_p_cut(dh::JuMP.Variable, q::Array{JuMP.Variable, 2}, dir::JuMP.Variable, q_sol::Float64, R::Array{Float64}, r_hat::Int, L::Float64)
    phi_hat = R[r_hat] * head_loss_hw_func(q_sol)
    phi_prime_hat = R[r_hat] * head_loss_hw_prime(q_sol)
    q_tilde = [compute_q_tilde(q_sol, R[r_hat], R[r]) for r in 1:length(R)]
    coeffs = [R[r] * head_loss_hw_prime(q_tilde[r]) for r in 1:length(R)]
    zero_indices = findall(r -> coeffs[r] < 1.0e-6, 1:length(R))
    setindex!(coeffs, zeros(size(zero_indices, 1)), zero_indices)

    expr = zero(AffExpr)
    expr += -dh / L + phi_prime_hat * sum(q[:, r_hat])
    expr += (phi_hat - phi_prime_hat * q_sol) * dir

    for r in setdiff(1:length(R), [r_hat])
        expr += coeffs[r] * sum(q[:, r])
    end

    return expr
end

function compute_q_n_cut(dh::JuMP.Variable, q::Array{JuMP.Variable, 2}, dir::JuMP.Variable, q_sol::Float64, R::Array{Float64}, r_hat::Int, L::Float64)
    phi_hat = R[r_hat] * head_loss_hw_func(q_sol)
    phi_prime_hat = R[r_hat] * head_loss_hw_prime(q_sol)
    q_tilde = [compute_q_tilde(q_sol, R[r_hat], R[r]) for r in 1:length(R)]
    coeffs = [R[r] * head_loss_hw_prime(q_tilde[r]) for r in 1:length(R)]
    zero_indices = findall(r -> coeffs[r] < 1.0e-6, 1:length(R))
    setindex!(coeffs, zeros(size(zero_indices, 1)), zero_indices)

    expr = zero(AffExpr)
    expr += -dh / L + phi_prime_hat * sum(q[:, r_hat])
    expr += (phi_hat - phi_prime_hat * q_sol) * (1 - dir)

    for r in setdiff(1:length(R), [r_hat])
        expr += coeffs[r] * sum(q[:, r])
    end

    return expr
end

function solution_is_integer(wm::GenericWaterModel, lp_solution::Array{Float64, 1}, n::Int = wm.cnw)
    dir_indices = linearindex.(wm.var[:nw][n][:dir][:])
    dir_sol = lp_solution[dir_indices]

    dir_zeros = filter(x -> isapprox(x, 0.0, atol = 0.01), dir_sol)
    dir_ones = filter(x -> isapprox(x, 1.0, atol = 0.01), dir_sol)
    num_dir_integer = length(dir_zeros) + length(dir_ones)

    if num_dir_integer < length(dir_indices)
        return false
    end

    for (a, connection) in wm.ref[:nw][n][:connection]
        xr_indices = linearindex.(wm.var[:nw][n][:xr][a])
        xr_sol = lp_solution[xr_indices]
        xr_zeros = filter(x -> isapprox(x, 0.0, atol = 0.01), xr_sol)
        xr_ones = filter(x -> isapprox(x, 1.0, atol = 0.01), xr_sol)
        num_xr_integer = length(xr_zeros) + length(xr_ones)

        if num_xr_integer < length(xr_indices)
            return false
        end
    end

    return true
end

function user_cut_callback_generator(wm::GenericWaterModel,
                                     params::Dict{String, Any},
                                     nlp_solver::MathProgBase.AbstractMathProgSolver,
                                     n::Int = wm.cnw)
    # Get indices associated with the MathProgBase model.
    arcs = collect(ids(wm, n, :connection))

    function user_cut_callback(cb::MathProgBase.MathProgCallbackData)
        lp_solution = MathProgBase.cbgetlpsolution(cb)
        dir_indices = linearindex.(wm.var[:nw][n][:dir][:])
        dir_sol = lp_solution[dir_indices]

        if solution_is_integer(wm, lp_solution, n)
            connection_ids = collect(ids(wm, n, :connection))
            resistance_indices = Dict{Int, Int}(a => 1 for a in connection_ids)

            for (a, connection) in wm.ref[:nw][n][:connection]
                xr_a = getvalue(wm.var[:nw][n][:xr][a])
                r = findfirst(i -> isapprox(xr_a[i], 1.0, atol = 0.01), 1:length(xr_a))
                resistance_indices[a] = r
            end

            # Update objective values.
            params["obj_last"] = params["obj_curr"]
            params["obj_curr"] = compute_objective(wm, resistance_indices, n)

            q, h = get_cvx_solution(wm, resistance_indices, nlp_solver)
            qlb, qub, hlb, hub = check_solution_bounds(wm, q, h, resistance_indices, n)
            solution_is_feasible = all([all(values(qlb)), all(values(qub)),
                                        all(values(hlb)), all(values(hub))])

            if solution_is_feasible
                for (a, connection) in wm.ref[:nw][n][:connection]
                    R_a = wm.ref[:nw][n][:resistance][a]
                    L_a = wm.ref[:nw][n][:connection][a]["length"]
                    dir = wm.var[:nw][n][:dir][a]
                    resistance_index = resistance_indices[a]

                    if q[a] >= 0.0
                        qp = wm.var[:nw][n][:qp][a]
                        dhp = wm.var[:nw][n][:dhp][a]
                        lhs = compute_q_p_cut(dhp, qp, dir, q[a], R_a, resistance_index, L_a)
                        @usercut(cb, lhs <= 0.0)
                    else
                        qn = wm.var[:nw][n][:qn][a]
                        dhn = wm.var[:nw][n][:dhn][a]
                        lhs = compute_q_n_cut(dhn, qn, dir, -q[a], R_a, resistance_index, L_a)
                        @usercut(cb, lhs <= 0.0)
                    end
                end

                params["obj_best"] = min(params["obj_best"], params["obj_curr"])
            end

            return
        end

        # Initialize the objective value.
        depth_satisfied = num_rounds_satisfied = true
        current_node = 0
        current_objective = 0.0

        for (a, connection) in wm.ref[:nw][n][:connection]
            L_a = wm.ref[:nw][n][:connection][a]["length"]
            C_a = wm.ref[:nw][n][:resistance_cost][a] .* L_a
            xr_ids = linearindex.(wm.var[:nw][n][:xr][a])
            current_objective += sum(lp_solution[xr_ids] .* C_a)
        end

        # Update objective values.
        params["obj_last"] = params["obj_curr"]
        params["obj_curr"] = current_objective

        # Conditions for adding outer approximations.
        obj_rel_change = (params["obj_curr"] - params["obj_last"]) / params["obj_last"]
        obj_improved = obj_rel_change >= params["K_oa"]

        if obj_improved
            for (relative_index, a) in enumerate(arcs)
                R_a = wm.ref[:nw][n][:resistance][a]
                L_a = wm.ref[:nw][n][:connection][a]["length"]
                dir = wm.var[:nw][n][:dir][a]

                if dir_sol[relative_index] >= 0.5
                    qp = wm.var[:nw][n][:qp][a]
                    qp_hat = lp_solution[linearindex.(qp)]
                    negative_indices = findall(signbit, qp_hat)
                    qp_hat[negative_indices] .= 0.0

                    dhp = wm.var[:nw][n][:dhp][a]
                    dhp_hat = lp_solution[linearindex(dhp)]
                    phi_max, r_hat = findmax([R_a[r] * sum(qp_hat[:, r].^(1.852)) for r in 1:length(R_a)])

                    if -dhp_hat / L_a + phi_max > params["epsilon"]
                        lhs = compute_q_p_cut(dhp, qp, dir, sum(qp_hat[:, r_hat]), R_a, r_hat, L_a)
                        @usercut(cb, lhs <= 0.0)
                    end
                else
                    qn = wm.var[:nw][n][:qn][a]
                    qn_hat = lp_solution[linearindex.(qn)]
                    negative_indices = findall(signbit, qn_hat)
                    qn_hat[negative_indices] .= 0.0

                    dhn = wm.var[:nw][n][:dhn][a]
                    dhn_hat = lp_solution[linearindex(dhn)]
                    phi_max, r_hat = findmax([R_a[r] * sum(qn_hat[:, r].^(1.852)) for r in 1:length(R_a)])

                    if -dhn_hat / L_a + phi_max > params["epsilon"]
                        lhs = compute_q_n_cut(dhn, qn, dir, sum(qn_hat[:, r_hat]), R_a, r_hat, L_a)
                        @usercut(cb, lhs <= 0.0)
                    end
                end
            end
        end
    end

    return user_cut_callback
end
