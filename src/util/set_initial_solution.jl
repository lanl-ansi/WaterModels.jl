function set_initial_solution(wm::GenericWaterModel, n_s::Int, R_id::Dict{Int, Int}, nlp_solver::MathProgBase.AbstractMathProgSolver)
    # Create the initial network.
    n = wm.cnw
    q, h = get_cvx_solution(wm, R_id, nlp_solver)

    for (i, junction) in wm.ref[:nw][n][:junctions]
        setvalue(wm.var[:nw][n][:h][i], h[i])
    end

    for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
        setvalue(wm.var[:nw][n][:h][i], h[i])
    end

    # Set resistances appropriately.
    for (a, connection) in wm.ref[:nw][n][:connection]
        dir = q[a] >= 0.0 ? 1 : 0
        i = parse(Int, connection["node1"])
        j = parse(Int, connection["node2"])

        setvalue.(wm.var[:nw][n][:qp][a], 0.0)
        setvalue.(wm.var[:nw][n][:qn][a], 0.0)
        setvalue.(wm.var[:nw][n][:xsn][a], 0)
        setvalue.(wm.var[:nw][n][:xsp][a], 0)
        setvalue.(wm.var[:nw][n][:xr][a], 0)

        setvalue(wm.var[:nw][n][:dir][a], dir)
        setvalue(wm.var[:nw][n][:xr][a][R_id[a]], 1)
        setvalue(wm.var[:nw][n][:dhp][a], dir * (h[i] - h[j]))
        setvalue(wm.var[:nw][n][:dhn][a], (1 - dir) * (h[j] - h[i]))

        for k in 1:n_s
            if q[a] >= 0.0
                qp_ar = wm.var[:nw][n][:qp][a][:, R_id[a]]
                qp_ark_lb = k > 1 ? getupperbound(qp_ar[k-1]) : 0.0
                qp_ark_ub = getupperbound(qp_ar[k])

                if q[a] >= qp_ark_lb && q[a] <= qp_ark_ub
                    setvalue(wm.var[:nw][n][:xsp][a][k, R_id[a]], 1)
                    setvalue(wm.var[:nw][n][:qp][a][k, R_id[a]], q[a])
                    break
                end
            else
                qn_ar = wm.var[:nw][n][:qn][a][:, R_id[a]]
                qn_ark_lb = k > 1 ? getupperbound(qn_ar[k-1]) : 0.0
                qn_ark_ub = getupperbound(qn_ar[k])

                if -q[a] >= qn_ark_lb && -q[a] <= qn_ark_ub
                    setvalue(wm.var[:nw][n][:xsn][a][k, R_id[a]], 1)
                    setvalue(wm.var[:nw][n][:qn][a][k, R_id[a]], -q[a])
                    break
                end
            end
        end
    end
end
