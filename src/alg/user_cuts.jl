export user_cut_callback_generator

import GLPK
import Random

function compute_q_tilde(q_hat::Float64, r_hat::Float64, r::Float64)
    if q_hat > 0.0
        phi_r_hat = r_hat * head_loss_hw_func(q_hat)
        phi_prime_r_hat = r_hat * head_loss_hw_prime(q_hat)
        b = phi_r_hat - phi_prime_r_hat
        return ((-(250 * b) / (213 * r))^(500 / 463))^(1 / 2)
    else
        return 0.0
    end
end

function compute_q_p_cut(dh::JuMP.Variable, q::Array{JuMP.Variable}, dir::JuMP.Variable, q_sol::Array{Float64}, R::Array{Float64}, r_hat::Int, L::Float64)
    phi_hat = R[r_hat] * head_loss_hw_func(q_sol[r_hat])
    phi_prime_hat = R[r_hat] * head_loss_hw_prime(q_sol[r_hat])
    q_tilde = [compute_q_tilde(q_sol[r_hat], R[r_hat], R[r]) for r in 1:length(R)]
    rneq = [r for r in 1:length(R) if r != r_hat]

    expr = zero(AffExpr)
    expr += -dh / L + phi_prime_hat * q[r_hat]
    expr += (phi_hat - phi_prime_hat * q_sol[r_hat]) * dir
    expr += sum(R[r] * head_loss_hw_prime(q_tilde[r]) * q[r] for r in rneq)
    return expr
end

function user_cut_callback_generator(wm::GenericWaterModel, params::Dict{String, Any}, n::Int = wm.cnw)
    # TODO: Not efficient... we need another method for storing resistances.
    R = calc_resistances_hw(wm, n)

    # Get indices associated with the MathProgBase model.
    arcs = collect(ids(wm, n, :connection))
    dir_indices = linearindex.(wm.var[:nw][n][:dir][:])
    xr_indices = vcat([linearindex.(wm.var[:nw][n][:xr][a][:]) for a in arcs]...)

    function user_cut_callback(cb::MathProgBase.MathProgCallbackData)
        lp_solution = MathProgBase.cbgetlpsolution(cb)
        dir_sol = lp_solution[dir_indices]
        xr_sol = lp_solution[xr_indices]

        # TODO: The method for finding these data should be solver-agnostic.
        current_node = GLPK.ios_curr_node(cb.tree)
        previous_node = GLPK.ios_prev_node(cb.tree, current_node)
        current_problem = GLPK.ios_get_prob(cb.tree)
        current_objective = GLPK.get_obj_val(current_problem)
        d = convert(Float64, GLPK.ios_node_level(cb.tree, current_node))

        # Update objective values.
        params["obj_last"] = params["obj_curr"]
        params["obj_curr"] = current_objective

        # Conditions for adding outer approximations.
        depth_satisfied = Random.rand() <= params["Beta_oa"] * 2^(-d)
        obj_change_rel = (params["obj_curr"] - params["obj_last"]) / params["obj_last"]
        obj_improved = obj_change_rel >= params["K_oa"]
       
        # Initialize the number of cut rounds added per node.
        if !haskey(params["n"], current_node)
            params["n"][current_node] = 0
        end

        # Check satisfaction of the number of rounds.
        num_rounds_satisfied = params["n"][current_node] <= params["M_oa"]

        if true #depth_satisfied && obj_improved && num_rounds_satisfied
            for (relative_index, a) in enumerate(arcs)
                L = wm.ref[:nw][n][:connection][a]["length"]
                dir = wm.var[:nw][n][:dir][a]

                if dir_sol[relative_index] >= 0.5
                    q_p_hat = lp_solution[linearindex.(wm.var[:nw][n][:qp][a][:])]
                    dhp_hat = lp_solution[linearindex(wm.var[:nw][n][:dhp][a])]
                    phi_max, r_hat = findmax([R[a][r] * q_p_hat[r]^(1.852) for r in 1:length(q_p_hat)])

                    if -dhp_hat / L + phi_max > params["epsilon"]
                        dh = wm.var[:nw][n][:dhp][a]
                        q = wm.var[:nw][n][:qp][a]
                        cut = compute_q_p_cut(dh, q, dir, q_p_hat, R[a], r_hat, L)
                        #@usercut(cb, cut <= 0.0)
                        #push!(params["oa_p"][a], (q_p_hat[r], r_hat))
                    end
                else
                    #q_n_hat = lp_solution[linearindex.(wm.var[:nw][n][:qn][a][:])]
                    #dhn_hat = lp_solution[linearindex(wm.var[:nw][n][:dhn][a])]
                    #phi_max, r_hat = findmax([R[a][r] * q_n_hat[r]^(1.852) for r in 1:length(q_n_hat)])

                    #if -dhn_hat / L + phi_max > params["epsilon"]
                    #    dh = wm.var[:nw][n][:dhn][a]
                    #    q = wm.var[:nw][n][:qn][a]
                    #    cut = compute_q_p_cut(dh, q, dir, q_n_hat, R[a], r_hat, L)
                    #    @usercut(cb, cut <= 0)
                    #    push!(params["oa_n"][a], (q_n_hat[r], r_hat))
                    #end
                end
            end

            params["n"][current_node] += 1
        end
    end

    return user_cut_callback
end
