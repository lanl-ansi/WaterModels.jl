function set_initial_solution(wm::GenericWaterModel, R_id::Dict{Int, Int}, nlp_solver::MathProgBase.AbstractMathProgSolver)
    # Create the initial network.
    n = wm.cnw
    q, h = get_cvx_solution(wm, R_id, nlp_solver)

    # Set resistances appropriately.
    for (a, connection) in wm.ref[:nw][n][:connection]
        setvalue.(wm.var[:nw][n][:qp][a], 0.0)
        setvalue(wm.var[:nw][n][:qp][a][R_id[a]], max(0.0, q[a]))

        setvalue.(wm.var[:nw][n][:qn][a], 0.0)
        setvalue(wm.var[:nw][n][:qn][a][R_id[a]], max(0.0, -q[a]))

        dir = q[a] > 0.0 ? 1 : 0
        setvalue(wm.var[:nw][n][:dir][a], dir)

        setvalue.(wm.var[:nw][n][:xr][a], 0.0)
        setvalue(wm.var[:nw][n][:xr][a][R_id[a]], 1.0)

        i = parse(Int, connection["node1"])
        j = parse(Int, connection["node2"])

        setvalue(wm.var[:nw][n][:dhp][a], dir * (h[i] - h[j]))
        setvalue(wm.var[:nw][n][:dhn][a], (1 - dir) * (h[j] - h[i]))
    end

    for (i, junction) in wm.ref[:nw][n][:junctions]
        setvalue(wm.var[:nw][n][:h][i], h[i])
    end

    for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
        setvalue(wm.var[:nw][n][:h][i], h[i])
    end
end
