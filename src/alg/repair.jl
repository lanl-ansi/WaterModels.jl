export repair_solution

import MathProgBase

"""
Implements an algorithm to repair a set of resistance choices whose network
analysis solution violates boundss on flow or potentials (Algorithm 3 in
Raghunathan (2013)). Returns a repaired set of resistance choices.
"""
function repair_solution(wm::GenericWaterModel, R_id::Dict{Int, Int}, params::Dict{String, Any}, nlp_solver::MathProgBase.AbstractMathProgSolver; feasibility_check = false)
    # Initialize parameters used within the while loop.
    n = wm.cnw
    num_iterations = 1
    progress = true
    repaired = false
    network = deepcopy(wm.data)

    # Initialize dictionaries used to store arc infeasibility results.
    connection_ids = collect(ids(wm, n, :connection))
    qp_sat_lb = Dict{Int, Bool}(a => true for a in connection_ids)
    qp_sat_ub = Dict{Int, Bool}(a => true for a in connection_ids)
    qn_sat_lb = Dict{Int, Bool}(a => true for a in connection_ids)
    qn_sat_ub = Dict{Int, Bool}(a => true for a in connection_ids)

    # Initialize dictionaries used to store node infeasibility results.
    junction_ids = collect(ids(wm, n, :junctions))
    reservoir_ids = collect(ids(wm, n, :reservoirs))
    node_ids = [junction_ids; reservoir_ids]
    h_sat_lb = Dict{Int, Bool}(i => true for i in node_ids)
    h_sat_ub = Dict{Int, Bool}(i => true for i in node_ids)

    while num_iterations <= params["max_repair_iters"] && !repaired && progress
        # Update resistances.
        for a in connection_ids
            selected_resistance = wm.ref[:nw][n][:resistance][a][R_id[a]]
            network["pipes"][string(a)]["resistance"] = selected_resistance
        end

        cvx = build_generic_model(network, CVXNLPWaterModel, WaterModels.post_cvx_hw)
        setsolver(cvx.model, nlp_solver)
        status = JuMP.solve(cvx.model)
        h = get_head_solution(cvx, n)

        for i in junction_ids
            h_sat_lb[i] = h[i] >= getlowerbound(wm.var[:nw][n][:h][i])
            h_sat_ub[i] = h[i] <= getupperbound(wm.var[:nw][n][:h][i])
        end

        for a in connection_ids
            qp = getvalue(cvx.var[:nw][n][:qp][a][1])
            qp_sat_lb[a] = qp >= getlowerbound(wm.var[:nw][n][:qp][a][R_id[a]])
            qp_sat_ub[a] = qp <= getupperbound(wm.var[:nw][n][:qp][a][R_id[a]])

            qn = getvalue(cvx.var[:nw][n][:qn][a][1])
            qn_sat_lb[a] = qn >= getlowerbound(wm.var[:nw][n][:qn][a][R_id[a]])
            qn_sat_ub[a] = qn <= getupperbound(wm.var[:nw][n][:qn][a][R_id[a]])
        end

        feasibility_checks = vcat([collect(values(h_sat_lb));
                                   collect(values(h_sat_ub));
                                   collect(values(qp_sat_lb));
                                   collect(values(qp_sat_ub));
                                   collect(values(qn_sat_lb));
                                   collect(values(qn_sat_ub))])

        if feasibility_check
            return all(feasibility_checks), R_id
        end

        if all(feasibility_checks)
            repaired = true
        else
            progress = false

            for a in connection_ids
                n_r = length(wm.ref[:nw][n][:resistance][a])

                if (!qp_sat_ub[a] || !qn_sat_ub[a]) && (R_id[a] < n_r)
                    R_id[a] += 1
                    progress = true
                end
            end

            if !progress
                for (a, connection) in wm.ref[:nw][n][:connection]
                    n_r = length(wm.ref[:nw][n][:resistance][a])

                    i = parse(Int, connection["node1"])
                    j = parse(Int, connection["node2"])

                    if (!h_sat_lb[j] && h_sat_lb[i]) || (!h_sat_lb[i] && h_sat_lb[j])
                       if R_id[a] < n_r
                           R_id[a] += 1
                           progress = true
                       end
                    end
                end
            end
        end

        # Initialize objective.
        objective = 0.0

        for (a, connection) in wm.ref[:nw][n][:connection]
            L_a = connection["length"]
            objective += wm.ref[:nw][n][:resistance_cost][a][R_id[a]] * L_a
        end

        if (objective >= params["obj_best"])
            progress = false
            repaired = false
        end

        num_iterations += 1
    end

    return repaired, R_id
end
