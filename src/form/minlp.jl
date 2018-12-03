# Defines MINLP implementations of water distribution models.

export MINLPWaterModel, StandardMINLPForm

abstract type AbstractMINLPForm <: AbstractWaterFormulation end
abstract type StandardMINLPForm <: AbstractMINLPForm end

"Default (nonconvex) MINLP model."
const MINLPWaterModel = GenericWaterModel{StandardMINLPForm}

"Default MINLP constructor."
MINLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMINLPForm; kwargs...)

function constraint_select_resistance(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: StandardMINLPForm
    if !haskey(wm.con[:nw][n], :select_resistance)
        wm.con[:nw][n][:select_resistance] = Dict{Int, ConstraintRef}()
    end

    con = @constraint(wm.model, sum(wm.var[:nw][n][:xr][a]) == 1)
    wm.con[:nw][n][:select_resistance] = con
end

function constraint_select_flow_term(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: StandardMINLPForm
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

    for r in keys(wm.var[:nw][n][:xr][a])
        q_p_r = wm.var[:nw][n][:qp][a][r]
        q_n_r = wm.var[:nw][n][:qn][a][r]
        x_dir = wm.var[:nw][n][:dir][a]
        x_r = wm.var[:nw][n][:xr][a][r]
        q_p_r_ub = getupperbound(q_p_r)
        q_n_r_ub = getupperbound(q_n_r)

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

function constraint_head_difference(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: StandardMINLPForm
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

function constraint_potential_loss(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: StandardMINLPForm
    if !haskey(wm.con[:nw][n], :potential_loss_1)
        wm.con[:nw][n][:potential_loss_1] = Dict{Int, Dict{Int, ConstraintRef}}()
        wm.con[:nw][n][:potential_loss_2] = Dict{Int, Dict{Int, ConstraintRef}}()
    end

    wm.con[:nw][n][:potential_loss_1][a] = Dict{Int, ConstraintRef}()
    wm.con[:nw][n][:potential_loss_2][a] = Dict{Int, ConstraintRef}()

    dhp = wm.var[:nw][n][:dhp][a]
    dhn = wm.var[:nw][n][:dhn][a]
    L = wm.ref[:nw][n][:connection][a]["length"]

    # TODO: Not efficient... we need another method for storing resistances.
    R = calc_resistances_hw(wm, n)

    for r in keys(wm.var[:nw][n][:xr][a])
        q_n_r = wm.var[:nw][n][:qn][a][r]
        q_p_r = wm.var[:nw][n][:qp][a][r]

        con_1 = @NLconstraint(wm.model, dhp / L >= R[a][r] * head_loss_hw(q_p_r))
        wm.con[:nw][n][:potential_loss_1][a][r] = con_1

        con_2 = @NLconstraint(wm.model, dhn / L >= R[a][r] * head_loss_hw(q_n_r))
        wm.con[:nw][n][:potential_loss_2][a][r] = con_2
    end
end

function constraint_potential_loss_slope(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: StandardMINLPForm
    if !haskey(wm.con[:nw][n], :potential_loss_slope_1)
        wm.con[:nw][n][:potential_loss_slope_1] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n][:potential_loss_slope_2] = Dict{Int, ConstraintRef}()
    end

    dhp = wm.var[:nw][n][:dhp][a]
    dhn = wm.var[:nw][n][:dhn][a]
    L = wm.ref[:nw][n][:connection][a]["length"]

    # TODO: Not efficient... we need another method for storing resistances.
    R = calc_resistances_hw(wm, n)

    qp_ubs = [getupperbound(wm.var[:nw][n][:qp][a][r]) for r in 1:length(R[a])]
    slopes_p = [R[a][r]*qp_ubs[r]^(1.852) / qp_ubs[r] for r in 1:length(R[a])]
    nan_indices = findall(isnan, slopes_p)
    setindex!(slopes_p, zeros(size(nan_indices, 1)), nan_indices)
    rhs_1 = AffExpr(wm.var[:nw][n][:qp][a], slopes_p, 0.0)
    con_1 = @constraint(wm.model, dhp / L <= rhs_1)
    wm.con[:nw][n][:potential_loss_slope_1] = con_1

    qn_ubs = [getupperbound(wm.var[:nw][n][:qn][a][r]) for r in 1:length(R[a])]
    slopes_n = [R[a][r]*qn_ubs[r]^(1.852) / qn_ubs[r] for r in 1:length(R[a])]
    nan_indices = findall(isnan, slopes_n)
    setindex!(slopes_n, zeros(size(nan_indices, 1)), nan_indices)
    rhs_2 = AffExpr(wm.var[:nw][n][:qn][a], slopes_n, 0.0)
    con_2 = @constraint(wm.model, dhn / L <= rhs_2)
    wm.con[:nw][n][:potential_loss_slope_2] = con_2
end

#function constraint_flow_direction(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: StandardMINLPForm
#    # Get source and target nodes corresponding to the edge.
#    i = parse(Int, wm.ref[:nw][n][:pipes][a]["node1"])
#    j = parse(Int, wm.ref[:nw][n][:pipes][a]["node2"])
#
#    # Collect variables needed for the constraint.
#    q = wm.var[:nw][n][:q][a]
#    h_i = wm.var[:nw][n][:h][i]
#    h_j = wm.var[:nw][n][:h][j]
#    dhp = wm.var[:nw][n][:dhp][a]
#    dhn = wm.var[:nw][n][:dhn][a]
#end

#"Non-convex Darcy-Weisbach constraint with unknown direction."
#function constraint_dw_unknown_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMINLPForm
#    # Collect variables and parameters needed for the constraint.
#    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)
#
#    # Add a nonlinear constraint for the head loss.
#    @NLconstraint(wm.model, h_i - h_j == lambda * q * abs(q))
#end
#
#"Non-convex Hazen-Williams constraint for flow with unknown direction."
#function constraint_hw_unknown_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMINLPForm
#    # Collect variables and parameters needed for the constraint.
#    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)
#
#    # Add a non-convex constraint for the head loss.
#    @NLconstraint(wm.model, h_i - h_j == lambda * q * (q^2)^0.426)
#end
#
