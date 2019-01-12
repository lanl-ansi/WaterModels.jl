function constraint_select_resistance(wm::GenericWaterModel, a::Int, n_n::Int)
    if !haskey(wm.con[:nw][n_n], :select_resistance)
        wm.con[:nw][n_n][:select_resistance] = Dict{Int, ConstraintRef}()
    end

    con = @constraint(wm.model, sum(wm.var[:nw][n_n][:xr][a]) == 1)
    wm.con[:nw][n_n][:select_resistance][a] = con
end

function constraint_select_segment(wm::GenericWaterModel, a::Int, n_n::Int)
    if !haskey(wm.con[:nw][n_n], :select_segment)
        wm.con[:nw][n_n][:select_segment] = Dict{Int, Dict{Int, ConstraintRef}}()
    end

    wm.con[:nw][n_n][:select_segment][a] = Dict{Int, ConstraintRef}()

    for r in 1:length(wm.ref[:nw][n_n][:resistance][a])
        x_r = wm.var[:nw][n_n][:xr][a][r]
        x_s_p = wm.var[:nw][n_n][:xsp][a][:, r]
        x_s_n = wm.var[:nw][n_n][:xsn][a][:, r]

        con = @constraint(wm.model, sum(x_s_p) + sum(x_s_n) == x_r)
        wm.con[:nw][n_n][:select_segment][a][r] = con
    end
end

function constraint_select_segmented_flow_term(wm::GenericWaterModel, a::Int, n_n::Int, n_s::Int)
    if !haskey(wm.con[:nw][n_n], :select_flow_term_1)
        wm.con[:nw][n_n][:select_flow_term_1] = Dict{Int, Dict{Int, Dict{Int, ConstraintRef}}}()
        wm.con[:nw][n_n][:select_flow_term_2] = Dict{Int, Dict{Int, Dict{Int, ConstraintRef}}}()
        wm.con[:nw][n_n][:select_flow_term_3] = Dict{Int, Dict{Int, Dict{Int, ConstraintRef}}}()
        wm.con[:nw][n_n][:select_flow_term_4] = Dict{Int, Dict{Int, Dict{Int, ConstraintRef}}}()
        wm.con[:nw][n_n][:select_flow_term_5] = Dict{Int, Dict{Int, Dict{Int, ConstraintRef}}}()
        wm.con[:nw][n_n][:select_flow_term_6] = Dict{Int, Dict{Int, Dict{Int, ConstraintRef}}}()
        wm.con[:nw][n_n][:select_flow_term_7] = Dict{Int, Dict{Int, Dict{Int, ConstraintRef}}}()
        wm.con[:nw][n_n][:select_flow_term_8] = Dict{Int, Dict{Int, Dict{Int, ConstraintRef}}}()
    end

    wm.con[:nw][n_n][:select_flow_term_1][a] = Dict{Int, Dict{Int, ConstraintRef}}()
    wm.con[:nw][n_n][:select_flow_term_2][a] = Dict{Int, Dict{Int, ConstraintRef}}()
    wm.con[:nw][n_n][:select_flow_term_3][a] = Dict{Int, Dict{Int, ConstraintRef}}()
    wm.con[:nw][n_n][:select_flow_term_4][a] = Dict{Int, Dict{Int, ConstraintRef}}()
    wm.con[:nw][n_n][:select_flow_term_5][a] = Dict{Int, Dict{Int, ConstraintRef}}()
    wm.con[:nw][n_n][:select_flow_term_6][a] = Dict{Int, Dict{Int, ConstraintRef}}()
    wm.con[:nw][n_n][:select_flow_term_7][a] = Dict{Int, Dict{Int, ConstraintRef}}()
    wm.con[:nw][n_n][:select_flow_term_8][a] = Dict{Int, Dict{Int, ConstraintRef}}()

    for k in 1:n_s
        wm.con[:nw][n_n][:select_flow_term_1][a][k] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n_n][:select_flow_term_2][a][k] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n_n][:select_flow_term_3][a][k] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n_n][:select_flow_term_4][a][k] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n_n][:select_flow_term_5][a][k] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n_n][:select_flow_term_6][a][k] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n_n][:select_flow_term_7][a][k] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n_n][:select_flow_term_8][a][k] = Dict{Int, ConstraintRef}()

            con_3 = @constraint(wm.model, q_p_akr <= q_p_ub_akr * x_dir)
            wm.con[:nw][n_n][:select_flow_term_3][a][k][r] = con_3

            con_4 = @constraint(wm.model, q_p_akr >= q_p_lb_akr * x_dir)
            wm.con[:nw][n_n][:select_flow_term_4][a][k][r] = con_4

            con_7 = @constraint(wm.model, q_n_akr <= q_n_ub_akr * (1 - x_dir))
            wm.con[:nw][n_n][:select_flow_term_7][a][k][r] = con_7

            con_8 = @constraint(wm.model, q_n_akr >= q_n_lb_akr * (1 - x_dir))
            wm.con[:nw][n_n][:select_flow_term_8][a][k][r] = con_8

        for r in 1:length(wm.ref[:nw][n_n][:resistance][a])
            q_n_lb_akr = q_p_lb_akr = 0.0
            q_n_akr = wm.var[:nw][n_n][:qn][a][k, r]
            q_n_ub_akr = getupperbound(q_n_akr)
            q_p_akr = wm.var[:nw][n_n][:qp][a][k, r]
            q_p_ub_akr = getupperbound(q_p_akr)

            if k > 1
                q_n_lb_akr = getupperbound(wm.var[:nw][n_n][:qn][a][k-1, r])
                q_p_lb_akr = getupperbound(wm.var[:nw][n_n][:qp][a][k-1, r])
            end

            x_dir = wm.var[:nw][n_n][:dir][a]
            x_r = wm.var[:nw][n_n][:xr][a][r]
            x_s_n = wm.var[:nw][n_n][:xsn][a][k, r]
            x_s_p = wm.var[:nw][n_n][:xsp][a][k, r]

            con_5 = @constraint(wm.model, q_p_akr <= q_p_ub_akr * x_s_p)
            wm.con[:nw][n_n][:select_flow_term_5][a][k][r] = con_1

            con_6 = @constraint(wm.model, q_p_akr >= q_p_lb_akr * x_s_p)
            wm.con[:nw][n_n][:select_flow_term_6][a][k][r] = con_2

            con_7 = @constraint(wm.model, q_n_akr <= q_n_ub_akr * x_s_n)
            wm.con[:nw][n_n][:select_flow_term_7][a][k][r] = con_5

            con_8 = @constraint(wm.model, q_n_akr >= q_n_lb_akr * x_s_n)
            wm.con[:nw][n_n][:select_flow_term_8][a][k][r] = con_6
        end
    end
end

function constraint_select_flow_term(wm::GenericWaterModel, a::Int, n_n::Int)
    if !haskey(wm.con[:nw][n_n], :select_flow_term_1)
        wm.con[:nw][n_n][:select_flow_term_1] = Dict{Int, Dict{Int, ConstraintRef}}()
        wm.con[:nw][n_n][:select_flow_term_2] = Dict{Int, Dict{Int, ConstraintRef}}()
        wm.con[:nw][n_n][:select_flow_term_3] = Dict{Int, Dict{Int, ConstraintRef}}()
        wm.con[:nw][n_n][:select_flow_term_4] = Dict{Int, Dict{Int, ConstraintRef}}()
    end

    wm.con[:nw][n_n][:select_flow_term_1][a] = Dict{Int, ConstraintRef}()
    wm.con[:nw][n_n][:select_flow_term_2][a] = Dict{Int, ConstraintRef}()
    wm.con[:nw][n_n][:select_flow_term_3][a] = Dict{Int, ConstraintRef}()
    wm.con[:nw][n_n][:select_flow_term_4][a] = Dict{Int, ConstraintRef}()

    for r in 1:length(wm.ref[:nw][n_n][:resistance][a])
        q_n_r = wm.var[:nw][n_n][:qn][a][r]
        q_n_r_ub = getupperbound(q_n_r)

        q_p_r = wm.var[:nw][n_n][:qp][a][r]
        q_p_r_ub = getupperbound(q_p_r)

        x_dir = wm.var[:nw][n_n][:dir][a]
        x_r = wm.var[:nw][n_n][:xr][a][r]

        con_1 = @constraint(wm.model, q_p_r <= q_p_r_ub * x_r)
        wm.con[:nw][n_n][:select_flow_term_1][a][r] = con_1

        con_2 = @constraint(wm.model, q_p_r <= q_p_r_ub * x_dir)
        wm.con[:nw][n_n][:select_flow_term_2][a][r] = con_2

        con_3 = @constraint(wm.model, q_n_r <= q_n_r_ub * x_r)
        wm.con[:nw][n_n][:select_flow_term_3][a][r] = con_3

        con_4 = @constraint(wm.model, q_n_r <= q_n_r_ub * (1 - x_dir))
        wm.con[:nw][n_n][:select_flow_term_4][a][r] = con_4
    end
end

function constraint_head_difference(wm::GenericWaterModel, a::Int, n_n::Int)
    if !haskey(wm.con[:nw][n_n], :head_difference_1)
        wm.con[:nw][n_n][:head_difference_1] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n_n][:head_difference_2] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n_n][:head_difference_3] = Dict{Int, ConstraintRef}()
    end

    x_dir = wm.var[:nw][n_n][:dir][a]
    dhn = wm.var[:nw][n_n][:dhn][a]
    dhp = wm.var[:nw][n_n][:dhp][a]
    dhn_max = getupperbound(dhn)
    dhp_max = getupperbound(dhp)

    i = parse(Int, wm.ref[:nw][n_n][:connection][a]["node1"])
    j = parse(Int, wm.ref[:nw][n_n][:connection][a]["node2"])
    h_i = wm.var[:nw][n_n][:h][i]
    h_j = wm.var[:nw][n_n][:h][j]

    con_1 = @constraint(wm.model, dhp <= dhp_max * x_dir)
    wm.con[:nw][n_n][:head_difference_1] = con_1

    con_2 = @constraint(wm.model, dhn <= dhn_max * (1 - x_dir))
    wm.con[:nw][n_n][:head_difference_2] = con_2

    con_3 = @constraint(wm.model, h_i - h_j == dhp - dhn)
    wm.con[:nw][n_n][:head_difference_3] = con_3
end

function constraint_potential_loss_slope_segmented(wm::GenericWaterModel, a::Int, n_n::Int, n_s::Int)
    if !haskey(wm.con[:nw][n_n], :potential_loss_slope_1)
        wm.con[:nw][n_n][:potential_loss_slope_1] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n_n][:potential_loss_slope_2] = Dict{Int, ConstraintRef}()
    end

    dhp = wm.var[:nw][n_n][:dhp][a]
    dhn = wm.var[:nw][n_n][:dhn][a]
    L = wm.ref[:nw][n_n][:connection][a]["length"]
    R_a = wm.ref[:nw][n_n][:resistance][a]

    qp_lbs = [getupperbound(wm.var[:nw][n_n][:qp][a][k, r]) for k in 1:n_s-1, r in 1:length(R_a)]
    qp_lbs = vcat(zeros((1, length(R_a))), qp_lbs)
    qp_ubs = [getupperbound(wm.var[:nw][n_n][:qp][a][k, r]) for k in 1:n_s, r in 1:length(R_a)]
    dqp = [qp_ubs[k, r] - qp_lbs[k, r] for k in 1:n_s, r in 1:length(R_a)]
    dfp = [R_a[r] * (qp_ubs[k, r]^(1.852) - qp_lbs[k, r]^(1.852)) for k in 1:n_s, r in 1:length(R_a)]
    qp_slopes = [dfp[k, r] / dqp[k, r] for k in 1:n_s, r in 1:length(R_a)]
    qp_constants = [R_a[r] * qp_lbs[k, r]^(1.852) for k in 1:n_s, r in 1:length(R_a)]
    qp_slopes = vcat(qp_slopes...)
    nan_indices = findall(isnan, qp_slopes)
    setindex!(qp_slopes, zeros(size(nan_indices, 1)), nan_indices)

    term_1 = AffExpr(vcat(wm.var[:nw][n_n][:qp][a]...), qp_slopes, 0.0)
    term_2 = AffExpr(vcat(wm.var[:nw][n_n][:xsp][a]...), vcat(qp_constants...), 0.0)
    
    con_1 = @constraint(wm.model, dhp / L <= term_1 + term_2)
    wm.con[:nw][n_n][:potential_loss_slope_1] = con_1

    qn_lbs = [getupperbound(wm.var[:nw][n_n][:qn][a][k, r]) for k in 1:n_s-1, r in 1:length(R_a)]
    qn_lbs = vcat(zeros((1, length(R_a))), qn_lbs)
    qn_ubs = [getupperbound(wm.var[:nw][n_n][:qn][a][k, r]) for k in 1:n_s, r in 1:length(R_a)]
    dqn = [qn_ubs[k, r] - qn_lbs[k, r] for k in 1:n_s, r in 1:length(R_a)]
    dfn = [R_a[r] * (qn_ubs[k, r]^(1.852) - qn_lbs[k, r]^(1.852)) for k in 1:n_s, r in 1:length(R_a)]
    qn_slopes = [dfn[k, r] / dqn[k, r] for k in 1:n_s, r in 1:length(R_a)]
    qn_constants = [R_a[r] * qn_lbs[k, r]^(1.852) for k in 1:n_s, r in 1:length(R_a)]
    qn_slopes = vcat(qn_slopes...)
    nan_indices = findall(isnan, qn_slopes)
    setindex!(qn_slopes, zeros(size(nan_indices, 1)), nan_indices)

    term_1 = AffExpr(vcat(wm.var[:nw][n_n][:qn][a]...), qn_slopes, 0.0)
    term_2 = AffExpr(vcat(wm.var[:nw][n_n][:xsn][a]...), vcat(qn_constants...), 0.0)
    
    con_2 = @constraint(wm.model, dhn / L <= term_1 + term_2)
    wm.con[:nw][n_n][:potential_loss_slope_2] = con_2
end

function constraint_potential_loss_slope(wm::GenericWaterModel, a::Int, n_n::Int)
    if !haskey(wm.con[:nw][n_n], :potential_loss_slope_1)
        wm.con[:nw][n_n][:potential_loss_slope_1] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n_n][:potential_loss_slope_2] = Dict{Int, ConstraintRef}()
    end

    dhp = wm.var[:nw][n_n][:dhp][a]
    dhn = wm.var[:nw][n_n][:dhn][a]
    L = wm.ref[:nw][n_n][:connection][a]["length"]
    R_a = wm.ref[:nw][n_n][:resistance][a]

    qp_ubs = [getupperbound(wm.var[:nw][n_n][:qp][a][r]) for r in 1:length(R_a)]
    slopes_p = [R_a[r]*qp_ubs[r]^(0.852) for r in 1:length(R_a)]
    rhs_1 = AffExpr(wm.var[:nw][n_n][:qp][a], slopes_p, 0.0)
    con_1 = @constraint(wm.model, dhp / L <= rhs_1)
    wm.con[:nw][n_n][:potential_loss_slope_1] = con_1

    qn_ubs = [getupperbound(wm.var[:nw][n_n][:qn][a][r]) for r in 1:length(R_a)]
    slopes_n = [R_a[r]*qn_ubs[r]^(0.852) for r in 1:length(R_a)]
    rhs_2 = AffExpr(wm.var[:nw][n_n][:qn][a], slopes_n, 0.0)
    con_2 = @constraint(wm.model, dhn / L <= rhs_2)
    wm.con[:nw][n_n][:potential_loss_slope_2] = con_2
end

"These problem forms use binary variables to specify flow direction."
AbstractDirectedForm = Union{AbstractMINLPForm, AbstractMILPRForm}

function variable_segment(wm::GenericWaterModel{T}, n_n::Int, n_s::Int) where T <: AbstractDirectedForm
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n_n, :connection))

    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    wm.var[:nw][n_n][:xsp] = Dict{Int, Array{Variable, 2}}()
    wm.var[:nw][n_n][:xsn] = Dict{Int, Array{Variable, 2}}()

    for (a, connection) in wm.ref[:nw][n_n][:connection]
        n_r = length(wm.ref[:nw][n_n][:resistance][a])

        for r in 1:n_r
            # Initialize variables associated with flow from i to j.
            wm.var[:nw][n_n][:xsp][a] = @variable(wm.model,
                                                  [k in 1:n_s, r in 1:n_r],
                                                  start = 0, category = :Bin,
                                                  basename = "xsp_$(n_n)_$(a)")

            wm.var[:nw][n_n][:xsn][a] = @variable(wm.model,
                                                  [k in 1:n_s, r in 1:n_r],
                                                  start = 0, category = :Bin,
                                                  basename = "xsn_$(n_n)_$(a)")

            setvalue(wm.var[:nw][n_n][:xsp][a][1, end], 1)
        end
    end
end

function variable_segmented_directed_flow(wm::GenericWaterModel{T}, n_n::Int, n_s::Int) where T <: AbstractDirectedForm
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n_n, :connection))

    # Compute sets of resistances.
    ub_n, ub_p = calc_directed_flow_upper_bounds(wm, n_n)

    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    wm.var[:nw][n_n][:qp] = Dict{Int, Array{Variable, 2}}()
    wm.var[:nw][n_n][:qn] = Dict{Int, Array{Variable, 2}}()

    for (a, connection) in wm.ref[:nw][n_n][:connection]
        n_r = length(wm.ref[:nw][n_n][:resistance][a])

        # Initialize variables associated with flow from i to j.
        wm.var[:nw][n_n][:qp][a] = @variable(wm.model, [k in 1:n_s, r in 1:n_r],
                                             lowerbound = 0.0,
                                             upperbound = k / n_s * ub_p[a][r],
                                             start = 0.0, category = :Cont,
                                             basename = "qp_$(n_n)_$(a)")

        # Initialize flow for the variable with least resistance.
        setvalue(wm.var[:nw][n_n][:qp][a][1, end], ub_p[a][1, end])

        # Initialize variables associated with flow from j to i.
        wm.var[:nw][n_n][:qn][a] = @variable(wm.model, [k in 1:n_s, r in 1:n_r],
                                             lowerbound = 0.0,
                                             upperbound = k/n_s * ub_n[a][r],
                                             start = 0.0, category = :Cont,
                                             basename = "qn_$(n_n)_$(a)")
    end
end

function variable_directed_flow(wm::GenericWaterModel{T}, n_n::Int) where T <: AbstractDirectedForm
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n_n, :connection))

    # Compute sets of resistances.
    ub_n, ub_p = calc_directed_flow_upper_bounds(wm, n_n)

    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    wm.var[:nw][n_n][:qp] = Dict{Int, Array{Variable, 1}}()
    wm.var[:nw][n_n][:qn] = Dict{Int, Array{Variable, 1}}()

    for (a, connection) in wm.ref[:nw][n_n][:connection]
        R_a = wm.ref[:nw][n_n][:resistance][a]

        # Initialize variables associated with flow from i to j.
        wm.var[:nw][n_n][:qp][a] = @variable(wm.model, [r in 1:length(R_a)],
                                             lowerbound = 0.0,
                                             upperbound = ub_p[a][r],
                                             start = 0.0, category = :Cont,
                                             basename = "qp_$(n_n)_$(a)")

        # Initialize flow for the variable with least resistance.
        setvalue(wm.var[:nw][n_n][:qp][a][1], ub_p[a][1])

        # Initialize variables associated with flow from j to i.
        wm.var[:nw][n_n][:qn][a] = @variable(wm.model, [r in 1:length(R_a)],
                                             lowerbound = 0.0,
                                             upperbound = ub_n[a][r],
                                             start = 0.0, category = :Cont,
                                             basename = "qn_$(n_n)_$(a)")
    end
end

"Constraint to ensure at least one direction is set to take flow away from a source."
function constraint_source_flow(wm::GenericWaterModel, i::Int, n::Int = wm.cnw)
    # Collect the required variables.
    connections = wm.ref[:nw][n][:connection]
    out_arcs = filter(a -> i == parse(Int, a.second["node1"]), connections)
    in_arcs = filter(a -> i == parse(Int, a.second["node2"]), connections)
    out_dirs = Array{JuMP.Variable}([wm.var[:nw][n][:dir][a] for a in keys(out_arcs)])
    in_dirs = Array{JuMP.Variable}([wm.var[:nw][n][:dir][a] for a in keys(in_arcs)])
    @constraint(wm.model, sum(out_dirs) - sum(in_dirs) >= 1 - length(in_dirs))
end

"Constraint to ensure at least one direction is set to take flow to a junction with demand."
function constraint_sink_flow(wm::GenericWaterModel, i::Int, n::Int = wm.cnw)
    # Collect the required variables.
    connections = wm.ref[:nw][n][:connection]
    out_arcs = filter(a -> i == parse(Int, a.second["node1"]), connections)
    in_arcs = filter(a -> i == parse(Int, a.second["node2"]), connections)
    out_dirs = Array{JuMP.Variable}([wm.var[:nw][n][:dir][a] for a in keys(out_arcs)])
    in_dirs = Array{JuMP.Variable}([wm.var[:nw][n][:dir][a] for a in keys(in_arcs)])
    @constraint(wm.model, sum(in_dirs) - sum(out_dirs) >= 1 - length(out_dirs))
end

#"Set new bounds for q given some specified direction of flow (-1 or 1)."
#function fix_flow_direction(q::JuMP.Variable, direction::Int)
#    # Fix the direction of the flow.
#    setlowerbound(q, direction == 1 ? 0.0 : getlowerbound(q))
#    setupperbound(q, direction == 1 ? getupperbound(q) : 0.0)
#end
#
#"Get variables commonly used in the construction of head loss constraints."
#function get_common_variables(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: AbstractWaterFormulation
#    # Get source and target nodes corresponding to the edge.
#    i = parse(Int, wm.ref[:nw][n][:pipes][a]["node1"])
#    j = parse(Int, wm.ref[:nw][n][:pipes][a]["node2"])
#
#    # Collect variables needed for the constraint.
#    q = wm.var[:nw][n][:q][a]
#    h_i = wm.var[:nw][n][:h][i]
#    h_j = wm.var[:nw][n][:h][j]
#
#    # Return the variables.
#    return q, h_i, h_j
#end
#
#"Get variables and constants used in the construction of Darcy-Weisbach constraints."
#function get_dw_requirements(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: AbstractWaterFormulation
#    q, h_i, h_j = get_common_variables(wm, a, n)
#    viscosity = wm.ref[:nw][n][:options]["viscosity"]
#    lambda = calc_friction_factor_dw(wm.ref[:nw][n][:pipes][a], viscosity)
#    return q, h_i, h_j, viscosity, lambda
#end
#
#"Get variables and constants used in the construction of Hazen-Williams constraints."
#function get_hw_requirements(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: AbstractWaterFormulation
#    q, h_i, h_j = get_common_variables(wm, a, n)
#    lambda = calc_friction_factor_hw(wm.ref[:nw][n][:pipes][a])
#    return q, h_i, h_j, lambda
#end
