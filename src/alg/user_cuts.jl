export user_cut_callback_generator, compute_q_tilde, compute_q_p_cut, compute_q_n_cut

#import CPLEX
import GLPK
import GLPKMathProgInterface
import Gurobi
import Random
import Roots

function compute_q_tilde(q_hat::Float64, r_hat::Float64, r::Float64)
    if q_hat > 0.0
        phi_r_hat = r_hat * head_loss_hw_func(q_hat)
        phi_prime_r_hat = r_hat * head_loss_hw_prime(q_hat)
        b = phi_r_hat - phi_prime_r_hat * q_hat
        return ((-(250.0 * b) / (213.0 * r))^(500.0 / 463.0))^(0.5)
    else
        return 0.0
    end
end

function compute_q_p_cut(dh::JuMP.Variable, q::Array{JuMP.Variable}, dir::JuMP.Variable, q_sol::Float64, R::Array{Float64}, r_hat::Int, L::Float64)
    phi_hat = R[r_hat] * head_loss_hw_func(q_sol)
    phi_prime_hat = R[r_hat] * head_loss_hw_prime(q_sol)
    q_tilde = [compute_q_tilde(q_sol, R[r_hat], R[r]) for r in 1:length(R)]
    rneq = setdiff(1:length(R), [r_hat])

    expr = zero(AffExpr)
    expr += -dh / L + phi_prime_hat * q[r_hat]
    expr += (phi_hat - phi_prime_hat * q_sol) * dir
    #expr += sum(R[r] * head_loss_hw_prime(q_tilde[r]) * q[r] for r in rneq)
    return expr
end

function compute_q_n_cut(dh::JuMP.Variable, q::Array{JuMP.Variable}, dir::JuMP.Variable, q_sol::Float64, R::Array{Float64}, r_hat::Int, L::Float64)
    phi_hat = R[r_hat] * head_loss_hw_func(q_sol)
    phi_prime_hat = R[r_hat] * head_loss_hw_prime(q_sol)
    q_tilde = [compute_q_tilde(q_sol, R[r_hat], R[r]) for r in 1:length(R)]
    rneq = setdiff(1:length(R), [r_hat])

    expr = zero(AffExpr)
    expr += -dh / L + phi_prime_hat * q[r_hat]
    expr += (phi_hat - phi_prime_hat * q_sol) * (1 - dir)
    #expr += sum(R[r] * head_loss_hw_prime(q_tilde[r]) * q[r] for r in rneq)
    return expr
end

function user_cut_callback_generator(wm::GenericWaterModel, params::Dict{String, Any}, n::Int = wm.cnw)
    # Get indices associated with the MathProgBase model.
    arcs = collect(ids(wm, n, :connection))
    dir_indices = linearindex.(wm.var[:nw][n][:dir][:])

    function user_cut_callback(cb::MathProgBase.MathProgCallbackData)
        lp_solution = MathProgBase.cbgetlpsolution(cb)
        dir_sol = lp_solution[dir_indices]

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

        if typeof(cb) == GLPKMathProgInterface.GLPKInterfaceMIP.GLPKCallbackData
            current_node = GLPK.ios_curr_node(cb.tree)
            current_problem = GLPK.ios_get_prob(cb.tree)
            d = convert(Float64, GLPK.ios_node_level(cb.tree, current_node))

            # Check satisfaction of node depth.
            depth_satisfied = Random.rand() <= params["Beta_oa"] * 2^(-d)

            # Initialize the number of cut rounds added per node.
            if !haskey(params["n"], current_node)
                params["n"][current_node] = 0
            end

            # Check satisfaction of the number of rounds.
            num_rounds_satisfied = params["n"][current_node] <= params["M_oa"]
        else
            # Initialize the number of cut rounds added per node.
            if !haskey(params["n"], current_node)
                params["n"][current_node] = 0
            end

            depth_satisfied = true
            num_rounds_satisfied = true
        end

        # Update objective values.
        params["obj_last"] = params["obj_curr"]
        params["obj_curr"] = current_objective

        # Conditions for adding outer approximations.
        obj_rel_change = (params["obj_curr"] - params["obj_last"]) / params["obj_last"]
        obj_improved = obj_rel_change >= params["K_oa"]

        if depth_satisfied && obj_improved && num_rounds_satisfied
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
                    phi_max, r_hat = findmax([R_a[r] * qp_hat[r]^(1.852) for r in 1:length(R_a)])

                    if -dhp_hat / L_a + phi_max > params["epsilon"]
                        lhs = compute_q_p_cut(dhp, qp, dir, qp_hat[r_hat], R_a, r_hat, L_a)
                        @usercut(cb, lhs <= 0.0)
                    end
                else
                    qn = wm.var[:nw][n][:qn][a]
                    qn_hat = lp_solution[linearindex.(qn)]
                    negative_indices = findall(signbit, qn_hat)
                    qn_hat[negative_indices] .= 0.0

                    dhn = wm.var[:nw][n][:dhn][a]
                    dhn_hat = lp_solution[linearindex(dhn)]
                    phi_max, r_hat = findmax([R_a[r] * qn_hat[r]^(1.852) for r in 1:length(R_a)])

                    if -dhn_hat / L_a + phi_max > params["epsilon"]
                        lhs = compute_q_n_cut(dhn, qn, dir, qn_hat[r_hat], R_a, r_hat, L_a)
                        @usercut(cb, lhs <= 0.0)
                    end
                end
            end

            params["n"][current_node] += 1
        end
    end

    return user_cut_callback
end
