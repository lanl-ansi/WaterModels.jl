function set_initial_solution(wm::GenericWaterModel, R_id::Dict{Int, Int}, nlp_solver::MathProgBase.AbstractMathProgSolver)
    # Create the initial network.
    n = wm.cnw
    network = deepcopy(wm.data)

    # Set resistances appropriately.
    for (a, connection) in wm.ref[:nw][n][:connection]
        selected_resistance = wm.ref[:nw][n][:resistance][a][R_id[a]]
        network["pipes"][string(a)]["resistance"] = selected_resistance
    end

    cvx = build_generic_model(network, CVXNLPWaterModel, WaterModels.post_cvx_hw)
    setsolver(cvx.model, nlp_solver)
    status = JuMP.solve(cvx.model)
    h = get_head_solution(cvx, n)

    # Set resistances appropriately.
    for (a, connection) in wm.ref[:nw][n][:connection]
        qp_sol = getvalue(cvx.var[:nw][n][:qp][a][1])
        setvalue.(wm.var[:nw][n][:qp][a], 0.0)
        setvalue(wm.var[:nw][n][:qp][a][R_id[a]], max(0.0, qp_sol))

        qn_sol = getvalue(cvx.var[:nw][n][:qn][a][1])
        setvalue.(wm.var[:nw][n][:qn][a], 0.0)
        setvalue(wm.var[:nw][n][:qn][a][R_id[a]], max(0.0, qn_sol))

        dir = (qp_sol - qn_sol) > 0.0 ? 1.0 : 0.0
        setvalue(wm.var[:nw][n][:dir][a], dir)

        setvalue.(wm.var[:nw][n][:xr][a], 0.0)
        setvalue(wm.var[:nw][n][:xr][a][R_id[a]], 1.0)

        h_i = h[parse(Int, connection["node1"])]
        h_j = h[parse(Int, connection["node2"])]

        setvalue(wm.var[:nw][n][:dhp][a], max(0.0, dir * (h_i - h_j)))
        setvalue(wm.var[:nw][n][:dhn][a], max(0.0, (dir - 1.0) * (h_i - h_j)))
    end

    for (i, junction) in wm.ref[:nw][n][:junctions]
        setvalue(wm.var[:nw][n][:h][i], h[i])
    end

    for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
        setvalue(wm.var[:nw][n][:h][i], h[i])
    end
end
