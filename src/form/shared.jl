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
    slopes_p = [R_a[r]*qp_ubs[r]^(1.852) / qp_ubs[r] for r in 1:length(R_a)]
    nan_indices = findall(isnan, slopes_p)
    setindex!(slopes_p, zeros(size(nan_indices, 1)), nan_indices)
    rhs_1 = AffExpr(wm.var[:nw][n][:qp][a], slopes_p, 0.0)
    con_1 = @constraint(wm.model, dhp / L <= rhs_1)
    wm.con[:nw][n][:potential_loss_slope_1] = con_1

    qn_ubs = [getupperbound(wm.var[:nw][n][:qn][a][r]) for r in 1:length(R_a)]
    slopes_n = [R_a[r]*qn_ubs[r]^(1.852) / qn_ubs[r] for r in 1:length(R_a)]
    nan_indices = findall(isnan, slopes_n)
    setindex!(slopes_n, zeros(size(nan_indices, 1)), nan_indices)
    rhs_2 = AffExpr(wm.var[:nw][n][:qn][a], slopes_n, 0.0)
    con_2 = @constraint(wm.model, dhn / L <= rhs_2)
    wm.con[:nw][n][:potential_loss_slope_2] = con_2
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
