#function constraint_select_resistance(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
#    if !haskey(wm.con[:nw][n], :select_resistance_sum)
#        wm.con[:nw][n][:select_resistance_sum] = Dict{Int, JuMP.ConstraintRef}()
#        wm.con[:nw][n][:select_resistance_n] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
#        wm.con[:nw][n][:select_resistance_p] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
#    end
#
#    con_sum = JuMP.@constraint(wm.model, sum(wm.var[:nw][n][:xr][a]) == 1)
#    wm.con[:nw][n][:select_resistance_sum][a] = con_sum
#
#    for r in 1:length(wm.ref[:nw][n][:resistance][a])
#        x_a_dir = wm.var[:nw][n][:dir][a]
#
#        q_n_a_r = wm.var[:nw][n][:qn][a][r]
#        q_n_a_r_ub = JuMP.upper_bound(q_n_a_r)
#        con_n = JuMP.@constraint(wm.model, q_n_a_r <= q_n_a_r_ub * (1 - x_dir))
#        wm.con[:nw][n][:select_flow_direction_n][a][r] = con_p
#
#        q_p_a_r = wm.var[:nw][n][:qp][a][r]
#        q_p_a_r_ub = JuMP.upper_bound(q_p_a_r)
#        con_p = JuMP.@constraint(wm.model, q_p_a_r <= q_p_a_r_ub * x_dir)
#        wm.con[:nw][n][:select_flow_direction_p][a][r] = con_p
#    end
#end
#
#function constraint_select_flow_direction(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
#    if !haskey(wm.con[:nw][n], :select_flow_direction_n)
#        wm.con[:nw][n][:select_flow_direction_n] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
#        wm.con[:nw][n][:select_flow_direction_p] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
#    end
#
#    wm.con[:nw][n][:select_flow_direction_n][a] = Dict{Int, JuMP.ConstraintRef}()
#    wm.con[:nw][n][:select_flow_direction_p][a] = Dict{Int, JuMP.ConstraintRef}()
#
#    for r in 1:length(wm.ref[:nw][n][:resistance][a])
#        x_a_res = wm.var[:nw][n][:xr][a][r]
#
#        q_n_a_r = wm.var[:nw][n][:qn][a][r]
#        q_n_a_r_ub = JuMP.upper_bound(q_n_a_r)
#        con_n = JuMP.@constraint(wm.model, q_n_a_r <= q_n_a_r_ub * x_a_res)
#        wm.con[:nw][n][:select_flow_direction_n][a][r] = con_n
#
#        q_p_a_r = wm.var[:nw][n][:qp][a][r]
#        q_p_a_r_ub = JuMP.upper_bound(q_p_a_r)
#        con_p = JuMP.@constraint(wm.model, q_p_a_r <= q_p_a_r_ub * x_a_res)
#        wm.con[:nw][n][:select_flow_direction_p][a][r] = con_p
#    end
#end

#function constraint_select_flow_term(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
#    if !haskey(wm.con[:nw][n], :select_flow_term_1)
#        wm.con[:nw][n][:select_flow_term_1] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
#        wm.con[:nw][n][:select_flow_term_2] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
#        wm.con[:nw][n][:select_flow_term_3] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
#        wm.con[:nw][n][:select_flow_term_4] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
#    end
#
#    wm.con[:nw][n][:select_flow_term_1][a] = Dict{Int, JuMP.ConstraintRef}()
#    wm.con[:nw][n][:select_flow_term_2][a] = Dict{Int, JuMP.ConstraintRef}()
#    wm.con[:nw][n][:select_flow_term_3][a] = Dict{Int, JuMP.ConstraintRef}()
#    wm.con[:nw][n][:select_flow_term_4][a] = Dict{Int, JuMP.ConstraintRef}()
#
#    for r in 1:length(wm.ref[:nw][n][:resistance][a])
#        q_n_r = wm.var[:nw][n][:qn][a][r]
#        q_n_r_ub = JuMP.upper_bound(q_n_r)
#
#        q_p_r = wm.var[:nw][n][:qp][a][r]
#        q_p_r_ub = JuMP.upper_bound(q_p_r)
#
#        x_dir = wm.var[:nw][n][:dir][a]
#        x_r = wm.var[:nw][n][:xr][a][r]
#
#        con_1 = JuMP.@constraint(wm.model, q_p_r <= q_p_r_ub * x_r)
#        wm.con[:nw][n][:select_flow_term_1][a][r] = con_1
#
#        con_2 = JuMP.@constraint(wm.model, q_p_r <= q_p_r_ub * x_dir)
#        wm.con[:nw][n][:select_flow_term_2][a][r] = con_2
#
#        con_3 = JuMP.@constraint(wm.model, q_n_r <= q_n_r_ub * x_r)
#        wm.con[:nw][n][:select_flow_term_3][a][r] = con_3
#
#        con_4 = JuMP.@constraint(wm.model, q_n_r <= q_n_r_ub * (1 - x_dir))
#        wm.con[:nw][n][:select_flow_term_4][a][r] = con_4
#    end
#end

#function constraint_head_difference(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
#    if !haskey(wm.con[:nw][n], :head_difference_1)
#        wm.con[:nw][n][:head_difference_1] = Dict{Int, JuMP.ConstraintRef}()
#        wm.con[:nw][n][:head_difference_2] = Dict{Int, JuMP.ConstraintRef}()
#        wm.con[:nw][n][:head_difference_3] = Dict{Int, JuMP.ConstraintRef}()
#    end
#
#    x_dir = wm.var[:nw][n][:dir][a]
#    dhn = wm.var[:nw][n][:dhn][a]
#    dhp = wm.var[:nw][n][:dhp][a]
#    dhn_max = JuMP.upper_bound(dhn)
#    dhp_max = JuMP.upper_bound(dhp)
#
#    i = wm.ref[:nw][n][:links][a]["node1"]
#    j = wm.ref[:nw][n][:links][a]["node2"]
#    h_i = wm.var[:nw][n][:h][i]
#    h_j = wm.var[:nw][n][:h][j]
#
#    con_1 = JuMP.@constraint(wm.model, dhp <= dhp_max * x_dir)
#    wm.con[:nw][n][:head_difference_1] = con_1
#
#    con_2 = JuMP.@constraint(wm.model, dhn <= dhn_max * (1 - x_dir))
#    wm.con[:nw][n][:head_difference_2] = con_2
#
#    con_3 = JuMP.@constraint(wm.model, h_i - h_j == dhp - dhn)
#    wm.con[:nw][n][:head_difference_3] = con_3
#end

#function constraint_potential_loss_slope(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
#    if !haskey(wm.con[:nw][n], :potential_loss_slope_1)
#        wm.con[:nw][n][:potential_loss_slope_1] = Dict{Int, JuMP.ConstraintRef}()
#        wm.con[:nw][n][:potential_loss_slope_2] = Dict{Int, JuMP.ConstraintRef}()
#    end
#
#    dhp = wm.var[:nw][n][:dhp][a]
#    dhn = wm.var[:nw][n][:dhn][a]
#    L = wm.ref[:nw][n][:links][a]["length"]
#    R_a = wm.ref[:nw][n][:resistance][a]
#
#    qp_ubs = [JuMP.upper_bound(wm.var[:nw][n][:qp][a][r]) for r in 1:length(R_a)]
#    slopes_p = [R_a[r]*qp_ubs[r]^(0.852) for r in 1:length(R_a)]
#    rhs_1 = sum(wm.var[:nw][n][:qp][a][r] * slopes_p[r] for r in 1:length(R_a))
#    con_1 = JuMP.@constraint(wm.model, dhp / L <= rhs_1)
#    wm.con[:nw][n][:potential_loss_slope_1] = con_1
#
#    qn_ubs = [JuMP.upper_bound(wm.var[:nw][n][:qn][a][r]) for r in 1:length(R_a)]
#    slopes_n = [R_a[r]*qn_ubs[r]^(0.852) for r in 1:length(R_a)]
#    rhs_2 = sum(wm.var[:nw][n][:qn][a][r] * slopes_n[r] for r in 1:length(R_a))
#    con_2 = JuMP.@constraint(wm.model, dhn / L <= rhs_2)
#    wm.con[:nw][n][:potential_loss_slope_2] = con_2
#end

"These problem forms use binary variables to specify flow direction."
AbstractDirectedForm = Union{AbstractMINLPForm, AbstractMILPRForm}

#function variable_directed_flow(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractDirectedForm
#    # Get indices for all network arcs.
#    arcs = collect(ids(wm, n, :links))
#
#    # Compute sets of resistances.
#    ub_n, ub_p = calc_directed_flow_upper_bounds(wm, n)
#
#    # Initialize directed flow variables. The variables qp correspond to flow
#    # from i to j, and the variables qn correspond to flow from j to i.
#    wm.var[:nw][n][:qp] = Dict{Int,Array{JuMP.VariableRef,1}}()
#    wm.var[:nw][n][:qn] = Dict{Int,Array{JuMP.VariableRef,1}}()
#
#    for (a, link) in wm.ref[:nw][n][:links]
#        R_a = wm.ref[:nw][n][:resistance][a]
#
#        # Initialize variables associated with flow from i to j.
#        wm.var[:nw][n][:qp][a] = JuMP.@variable(wm.model, [r in 1:length(R_a)],
#                                                lower_bound = 0.0,
#                                                upper_bound = ub_p[a][r],
#                                                start = 0.0,
#                                                base_name = "q_p[$(n)][$(a)]")
#
#        # Initialize variables associated with flow from j to i.
#        wm.var[:nw][n][:qn][a] = JuMP.@variable(wm.model, [r in 1:length(R_a)],
#                                                lower_bound = 0.0,
#                                                upper_bound = ub_n[a][r],
#                                                start = 0.0,
#                                                base_name = "q_n[$(n)][$(a)]")
#
#        # Initialize flow for the variable with least resistance.
#        JuMP.set_start_value(wm.var[:nw][n][:qp][a][1], ub_p[a][1])
#    end
#end

#"Constraint to ensure at least one direction is set to take flow away from a source."
#function constraint_source_flow(wm::GenericWaterModel, i::Int, n::Int=wm.cnw)
#    # Collect the required variables.
#    links = wm.ref[:nw][n][:links]
#    out_arcs = filter(a -> i == a.second["node1"], links)
#    in_arcs = filter(a -> i == a.second["node2"], links)
#    out_dirs = Array{JuMP.VariableRef}([wm.var[:nw][n][:dir][a] for a in keys(out_arcs)])
#    in_dirs = Array{JuMP.VariableRef}([wm.var[:nw][n][:dir][a] for a in keys(in_arcs)])
#    JuMP.@constraint(wm.model, sum(out_dirs) - sum(in_dirs) >= 1 - length(in_dirs))
#end
#
#"Constraint to ensure at least one direction is set to take flow to a junction with demand."
#function constraint_sink_flow(wm::GenericWaterModel, i::Int, n::Int=wm.cnw)
#    # Collect the required variables.
#    links = wm.ref[:nw][n][:links]
#    out_arcs = filter(a -> i == a.second["node1"], links)
#    in_arcs = filter(a -> i == a.second["node2"], links)
#    out_dirs = Array{JuMP.VariableRef}([wm.var[:nw][n][:dir][a] for a in keys(out_arcs)])
#    in_dirs = Array{JuMP.VariableRef}([wm.var[:nw][n][:dir][a] for a in keys(in_arcs)])
#    JuMP.@constraint(wm.model, sum(in_dirs) - sum(out_dirs) >= 1 - length(out_dirs))
#end

#"Set new bounds for q given some specified direction of flow (-1 or 1)."
#function fix_flow_direction(q::JuMP.Variable, direction::Int)
#    # Fix the direction of the flow.
#    setlower_bound(q, direction == 1 ? 0.0 : lower_bound(q))
#    setupper_bound(q, direction == 1 ? upper_bound(q) : 0.0)
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
