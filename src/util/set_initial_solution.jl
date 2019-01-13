function set_initial_solution(wm::GenericWaterModel, n_s::Int, R_id::Dict{Int, Int}, nlp_solver::MathProgBase.AbstractMathProgSolver)
    # Create the initial network.
    n = wm.cnw
    q, h = get_cvx_solution(wm, R_id, nlp_solver)

    # Set resistances appropriately.
    for (a, connection) in wm.ref[:nw][n][:connection]
        setvalue.(wm.var[:nw][n][:qp][a], 0.0)
        setvalue.(wm.var[:nw][n][:qn][a], 0.0)
        setvalue.(wm.var[:nw][n][:xsp][a], 0.0)
        setvalue.(wm.var[:nw][n][:xsn][a], 0.0)

        for k in 1:n_s
            if q[a] >= 0.0
                q_lb = 0.0

                if k > 1
                    q_lb = getupperbound(wm.var[:nw][n][:qp][a][k-1, R_id[a]])
                end

                q_ub = getupperbound(wm.var[:nw][n][:qp][a][k, R_id[a]])

                if q[a] >= q_lb && q[a] <= q_ub
                    setvalue(wm.var[:nw][n][:qp][a][k, R_id[a]], q[a])
                    setvalue(wm.var[:nw][n][:xsp][a][k, R_id[a]], 1)
                end
            else
                q_lb = 0.0

                if k > 1
                    q_lb = getupperbound(wm.var[:nw][n][:qn][a][k-1, R_id[a]])
                end

                q_ub = getupperbound(wm.var[:nw][n][:qn][a][k, R_id[a]])

                if -q[a] >= q_lb && -q[a] <= q_ub
                    setvalue(wm.var[:nw][n][:qn][a][k, R_id[a]], -q[a])
                    setvalue(wm.var[:nw][n][:xsn][a][k, R_id[a]], 1)
                end
            end
        end

        dir = q[a] > 0.0 ? 1 : 0
        setvalue(wm.var[:nw][n][:dir][a], dir)

        setvalue.(wm.var[:nw][n][:xr][a], 0)
        setvalue(wm.var[:nw][n][:xr][a][R_id[a]], 1)

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
