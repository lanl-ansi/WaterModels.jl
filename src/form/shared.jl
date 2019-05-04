function constraint_select_resistance(wm::GenericWaterModel, a::Int, n::Int = wm.cnw)
    if !haskey(wm.con[:nw][n], :select_resistance)
        wm.con[:nw][n][:select_resistance] = Dict{Int, ConstraintRef}()
    end

    con = @constraint(wm.model, sum(wm.var[:nw][n][:xr][a]) == 1)
    wm.con[:nw][n][:select_resistance][a] = con
end

function constraint_select_flow_term(wm::GenericWaterModel, a::Int, n::Int = wm.cnw)
    if !haskey(wm.con[:nw][n], :select_flow_term_1)
        wm.con[:nw][n][:select_flow_term_1] = Dict{Int, Dict{Int, ConstraintRef}}()
        wm.con[:nw][n][:select_flow_term_2] = Dict{Int, Dict{Int, ConstraintRef}}()
        wm.con[:nw][n][:select_flow_term_3] = Dict{Int, Dict{Int, ConstraintRef}}()
        wm.con[:nw][n][:select_flow_term_4] = Dict{Int, Dict{Int, ConstraintRef}}()
    end

    wm.con[:nw][n][:select_flow_term_1][a] = Dict{Int, ConstraintRef}()
    wm.con[:nw][n][:select_flow_term_2][a] = Dict{Int, ConstraintRef}()
    wm.con[:nw][n][:select_flow_term_3][a] = Dict{Int, ConstraintRef}()
    wm.con[:nw][n][:select_flow_term_4][a] = Dict{Int, ConstraintRef}()

    for r in 1:length(wm.ref[:nw][n][:resistance][a])
        q_n_r = wm.var[:nw][n][:qn][a][r]
        q_n_r_ub = getupperbound(q_n_r)

        q_p_r = wm.var[:nw][n][:qp][a][r]
        q_p_r_ub = getupperbound(q_p_r)

        x_dir = wm.var[:nw][n][:dir][a]
        x_r = wm.var[:nw][n][:xr][a][r]

        con_1 = @constraint(wm.model, q_p_r <= q_p_r_ub * x_r)
        wm.con[:nw][n][:select_flow_term_1][a][r] = con_1

        con_2 = @constraint(wm.model, q_p_r <= q_p_r_ub * x_dir)
        wm.con[:nw][n][:select_flow_term_2][a][r] = con_2

        con_3 = @constraint(wm.model, q_n_r <= q_n_r_ub * x_r)
        wm.con[:nw][n][:select_flow_term_3][a][r] = con_3

        con_4 = @constraint(wm.model, q_n_r <= q_n_r_ub * (1 - x_dir))
        wm.con[:nw][n][:select_flow_term_4][a][r] = con_4
    end
end

function constraint_head_difference(wm::GenericWaterModel, a::Int, n::Int = wm.cnw)
    if !haskey(wm.con[:nw][n], :head_difference_1)
        wm.con[:nw][n][:head_difference_1] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n][:head_difference_2] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n][:head_difference_3] = Dict{Int, ConstraintRef}()
    end

    x_dir = wm.var[:nw][n][:dir][a]
    dhn = wm.var[:nw][n][:dhn][a]
    dhp = wm.var[:nw][n][:dhp][a]
    dhn_max = getupperbound(dhn)
    dhp_max = getupperbound(dhp)

    i = parse(Int, wm.ref[:nw][n][:connection][a]["node1"])
    j = parse(Int, wm.ref[:nw][n][:connection][a]["node2"])
    h_i = wm.var[:nw][n][:h][i]
    h_j = wm.var[:nw][n][:h][j]

    con_1 = @constraint(wm.model, dhp <= dhp_max * x_dir)
    wm.con[:nw][n][:head_difference_1] = con_1

    con_2 = @constraint(wm.model, dhn <= dhn_max * (1 - x_dir))
    wm.con[:nw][n][:head_difference_2] = con_2

    con_3 = @constraint(wm.model, h_i - h_j == dhp - dhn)
    wm.con[:nw][n][:head_difference_3] = con_3
end

function constraint_potential_loss_slope(wm::GenericWaterModel, a::Int, n::Int = wm.cnw)
    if !haskey(wm.con[:nw][n], :potential_loss_slope_1)
        wm.con[:nw][n][:potential_loss_slope_1] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n][:potential_loss_slope_2] = Dict{Int, ConstraintRef}()
    end

    dhp = wm.var[:nw][n][:dhp][a]
    dhn = wm.var[:nw][n][:dhn][a]
    L = wm.ref[:nw][n][:connection][a]["length"]
    R_a = wm.ref[:nw][n][:resistance][a]

    qp_ubs = [getupperbound(wm.var[:nw][n][:qp][a][r]) for r in 1:length(R_a)]
    slopes_p = [R_a[r]*qp_ubs[r]^(0.852) for r in 1:length(R_a)]
    rhs_1 = AffExpr(wm.var[:nw][n][:qp][a], slopes_p, 0.0)
    con_1 = @constraint(wm.model, dhp / L <= rhs_1)
    wm.con[:nw][n][:potential_loss_slope_1] = con_1

    qn_ubs = [getupperbound(wm.var[:nw][n][:qn][a][r]) for r in 1:length(R_a)]
    slopes_n = [R_a[r]*qn_ubs[r]^(0.852) for r in 1:length(R_a)]
    rhs_2 = AffExpr(wm.var[:nw][n][:qn][a], slopes_n, 0.0)
    con_2 = @constraint(wm.model, dhn / L <= rhs_2)
    wm.con[:nw][n][:potential_loss_slope_2] = con_2
end

"These problem forms use binary variables to specify flow direction."
AbstractDirectedForm = Union{AbstractMICPForm, AbstractMILPRForm}

function variable_directed_flow(wm::GenericWaterModel{T}, n::Int = wm.cnw) where T <: AbstractDirectedForm
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n, :connection))

    # Compute sets of resistances.
    ub_n, ub_p = calc_directed_flow_upper_bounds(wm, n)

    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    wm.var[:nw][n][:qp] = Dict{Int, Array{Variable, 1}}()
    wm.var[:nw][n][:qn] = Dict{Int, Array{Variable, 1}}()

    for (a, connection) in wm.ref[:nw][n][:connection]
        R_a = wm.ref[:nw][n][:resistance][a]

        # Initialize variables associated with flow from i to j.
        wm.var[:nw][n][:qp][a] = @variable(wm.model, [r in 1:length(R_a)],
                                             lowerbound = 0.0,
                                             upperbound = ub_p[a][r],
                                             start = 0.0, category = :Cont,
                                             basename = "qp_$(n)_$(a)")

        # Initialize variables associated with flow from j to i.
        wm.var[:nw][n][:qn][a] = @variable(wm.model, [r in 1:length(R_a)],
                                             lowerbound = 0.0,
                                             upperbound = ub_n[a][r],
                                             start = 0.0, category = :Cont,
                                             basename = "qn_$(n)_$(a)")

        # Initialize flow for the variable with least resistance.
        setvalue(wm.var[:nw][n][:qp][a][1], ub_p[a][1])
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
