########################################################################
# This file defines commonly-used constraints for water systems models.
########################################################################

function constraint_undirected_flow_conservation(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Create the constraint dictionary if necessary.
    if !haskey(wm.con[:nw][n], :flow_conservation)
        wm.con[:nw][n][:flow_conservation] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Initialize the flow sum expression.
    flow_sum = JuMP.AffExpr(0.0)

    # Add all incoming flow to node i.
    for (a, link) in filter(is_in_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, wm.var[:nw][n][:q][a])
    end

    # Subtract all outgoing flow from node i.
    for (a, link) in filter(is_out_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, -wm.var[:nw][n][:q][a])
    end

    # Add the flow conservation constraint.
    demand = wm.ref[:nw][n][:junctions][i]["demand"]
    con = JuMP.@constraint(wm.model, flow_sum == demand)
    wm.con[:nw][n][:flow_conservation][i] = con
end

function constraint_undirected_flow_conservation_ne(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Create the constraint dictionary if necessary.
    if !haskey(wm.con[:nw][n], :flow_conservation_ne)
        wm.con[:nw][n][:flow_conservation_ne] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Initialize the flow sum expression.
    flow_sum = JuMP.AffExpr(0.0)

    # Add all incoming flow to node i.
    for (a, conn) in filter(is_in_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, sum(wm.var[:nw][n][:q_ne][a]))
    end

    # Subtract all outgoing flow from node i.
    for (a, conn) in filter(is_out_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, -sum(wm.var[:nw][n][:q_ne][a]))
    end

    # Add the flow conservation constraint.
    demand = wm.ref[:nw][n][:junctions][i]["demand"]
    con = JuMP.@constraint(wm.model, flow_sum == demand)
    wm.con[:nw][n][:flow_conservation_ne][i] = con
end

function constraint_directed_flow_conservation(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Create the constraint dictionary if necessary.
    if !haskey(wm.con[:nw][n], :flow_conservation)
        wm.con[:nw][n][:flow_conservation] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Initialize the flow sum expression.
    flow_sum = JuMP.AffExpr(0.0)

    # Add all incoming flow to node i.
    for (a, conn) in filter(is_in_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, -wm.var[:nw][n][:qn][a])
        JuMP.add_to_expression!(flow_sum, wm.var[:nw][n][:qp][a])
    end

    # Subtract all outgoing flow from node i.
    for (a, conn) in filter(is_out_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, wm.var[:nw][n][:qn][a])
        JuMP.add_to_expression!(flow_sum, -wm.var[:nw][n][:qp][a])
    end

    # Add the flow conservation constraint.
    demand = wm.ref[:nw][n][:junctions][i]["demand"]
    con = JuMP.@constraint(wm.model, flow_sum == demand)
    wm.con[:nw][n][:flow_conservation][i] = con
end

function constraint_directed_flow_conservation_ne(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Create the constraint dictionary if necessary.
    if !haskey(wm.con[:nw][n], :flow_conservation_ne)
        wm.con[:nw][n][:flow_conservation_ne] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Initialize the flow sum expression.
    flow_sum = JuMP.AffExpr(0.0)

    # Add all incoming flow to node i.
    for (a, conn) in filter(is_in_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, -sum(wm.var[:nw][n][:qn_ne][a]))
        JuMP.add_to_expression!(flow_sum, sum(wm.var[:nw][n][:qp_ne][a]))
    end

    # Subtract all outgoing flow from node i.
    for (a, conn) in filter(is_out_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, sum(wm.var[:nw][n][:qn_ne][a]))
        JuMP.add_to_expression!(flow_sum, -sum(wm.var[:nw][n][:qp_ne][a]))
    end

    # Add the flow conservation constraint.
    demand = wm.ref[:nw][n][:junctions][i]["demand"]
    con = JuMP.@constraint(wm.model, flow_sum == demand)
    wm.con[:nw][n][:flow_conservation_ne][i] = con
end

function constraint_directed_resistance_selection_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    if !haskey(wm.con[:nw][n], :resistance_selection_sum)
        wm.con[:nw][n][:resistance_selection_sum] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:resistance_selection_n] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        wm.con[:nw][n][:resistance_selection_p] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    wm.con[:nw][n][:resistance_selection_n][a] = Dict{Int, JuMP.ConstraintRef}()
    wm.con[:nw][n][:resistance_selection_p][a] = Dict{Int, JuMP.ConstraintRef}()

    con_sum = JuMP.@constraint(wm.model, sum(wm.var[:nw][n][:x_res][a]) == 1.0)
    wm.con[:nw][n][:resistance_selection_sum][a] = con_sum

    for r in 1:length(wm.ref[:nw][n][:resistance][a])
        x_res = wm.var[:nw][n][:x_res][a][r]

        qn_ne = wm.var[:nw][n][:qn_ne][a][r]
        qn_ne_ub = JuMP.upper_bound(qn_ne)
        con_n = JuMP.@constraint(wm.model, qn_ne - qn_ne_ub * x_res <= 0.0)
        wm.con[:nw][n][:resistance_selection_n][a][r] = con_n

        qp_ne = wm.var[:nw][n][:qp_ne][a][r]
        qp_ne_ub = JuMP.upper_bound(qp_ne)
        con_p = JuMP.@constraint(wm.model, qp_ne - qp_ne_ub * x_res <= 0.0)
        wm.con[:nw][n][:resistance_selection_p][a][r] = con_p
    end
end

function constraint_undirected_resistance_selection_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    if !haskey(wm.con[:nw][n], :resistance_selection_sum)
        wm.con[:nw][n][:resistance_selection_sum] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:resistance_selection_lb] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        wm.con[:nw][n][:resistance_selection_ub] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    con_sum = JuMP.@constraint(wm.model, sum(wm.var[:nw][n][:x_res][a]) == 1.0)
    wm.con[:nw][n][:resistance_selection_sum][a] = con_sum

    wm.con[:nw][n][:resistance_selection_lb][a] = Dict{Int, JuMP.ConstraintRef}()
    wm.con[:nw][n][:resistance_selection_ub][a] = Dict{Int, JuMP.ConstraintRef}()

    for r in 1:length(wm.ref[:nw][n][:resistance][a])
        x_res = wm.var[:nw][n][:x_res][a][r]

        q_ne = wm.var[:nw][n][:q_ne][a][r]
        q_ne_lb = JuMP.lower_bound(q_ne)
        con_lb = JuMP.@constraint(wm.model, q_ne - q_ne_lb * x_res >= 0.0)
        wm.con[:nw][n][:resistance_selection_lb][a][r] = con_lb

        q_ne = wm.var[:nw][n][:q_ne][a][r]
        q_ne_ub = JuMP.upper_bound(q_ne)
        con_ub = JuMP.@constraint(wm.model, q_ne - q_ne_ub * x_res <= 0.0)
        wm.con[:nw][n][:resistance_selection_ub][a][r] = con_ub
    end
end

function constraint_flow_direction_selection(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :flow_direction_selection_n)
        wm.con[:nw][n][:flow_direction_selection_n] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:flow_direction_selection_p] = Dict{Int, JuMP.ConstraintRef}()
    end

    x_dir = wm.var[:nw][n][:x_dir][a]

    qn = wm.var[:nw][n][:qn][a]
    qn_ub = JuMP.upper_bound(qn)
    con_n = JuMP.@constraint(wm.model, qn - qn_ub * (1.0 - x_dir) <= 0.0)
    wm.con[:nw][n][:flow_direction_selection_n][a] = con_n

    qp = wm.var[:nw][n][:qp][a]
    qp_ub = JuMP.upper_bound(qp)
    con_p = JuMP.@constraint(wm.model, qp - qp_ub * x_dir <= 0.0)
    wm.con[:nw][n][:flow_direction_selection_p][a] = con_p
end

function constraint_flow_direction_selection_ne(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :flow_direction_selection_ne_n)
        wm.con[:nw][n][:flow_direction_selection_ne_n] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        wm.con[:nw][n][:flow_direction_selection_ne_p] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    wm.con[:nw][n][:flow_direction_selection_ne_n][a] = Dict{Int, JuMP.ConstraintRef}()
    wm.con[:nw][n][:flow_direction_selection_ne_p][a] = Dict{Int, JuMP.ConstraintRef}()

    for r in 1:length(wm.ref[:nw][n][:resistance][a])
        x_dir = wm.var[:nw][n][:x_dir][a]

        qn_ne = wm.var[:nw][n][:qn_ne][a][r]
        qn_ne_ub = JuMP.upper_bound(qn_ne)
        con_n = JuMP.@constraint(wm.model, qn_ne - qn_ne_ub * (1.0 - x_dir) <= 0.0)
        wm.con[:nw][n][:flow_direction_selection_ne_n][a][r] = con_n

        qp_ne = wm.var[:nw][n][:qp_ne][a][r]
        qp_ne_ub = JuMP.upper_bound(qp_ne)
        con_p = JuMP.@constraint(wm.model, qp_ne - qp_ne_ub * x_dir <= 0.0)
        wm.con[:nw][n][:flow_direction_selection_ne_p][a][r] = con_p
    end
end

function constraint_directed_head_difference(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :head_difference_1)
        wm.con[:nw][n][:head_difference_1] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:head_difference_2] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:head_difference_3] = Dict{Int, JuMP.ConstraintRef}()
    end

    i = wm.ref[:nw][n][:links][a]["f_id"]

    if i in collect(ids(wm, n, :reservoirs))
        h_i = wm.ref[:nw][n][:reservoirs][i]["head"]
    else
        h_i = wm.var[:nw][n][:h][i]
    end

    j = wm.ref[:nw][n][:links][a]["t_id"]

    if j in collect(ids(wm, n, :reservoirs))
        h_j = wm.ref[:nw][n][:reservoirs][j]["head"]
    else
        h_j = wm.var[:nw][n][:h][j]
    end

    x_dir = wm.var[:nw][n][:x_dir][a]

    dhn = wm.var[:nw][n][:dhn][a]
    dhn_ub = JuMP.upper_bound(dhn)
    con_1 = JuMP.@constraint(wm.model, dhn - dhn_ub * (1.0 - x_dir) <= 0.0)
    wm.con[:nw][n][:head_difference_1][a] = con_1

    dhp = wm.var[:nw][n][:dhp][a]
    dhp_ub = JuMP.upper_bound(dhp)
    con_2 = JuMP.@constraint(wm.model, dhp - dhp_ub * x_dir <= 0.0)
    wm.con[:nw][n][:head_difference_2][a] = con_2

    con_3 = JuMP.@constraint(wm.model, (h_i - h_j) - (dhp - dhn) == 0.0)
    wm.con[:nw][n][:head_difference_3][a] = con_3
end

function constraint_directed_potential_loss_ub_ne(wm::GenericWaterModel, a::Int, n::Int=wm.cnw; alpha::Float64=1.852)
    if !haskey(wm.con[:nw][n], :directed_potential_loss_ub_ne_n)
        wm.con[:nw][n][:directed_potential_loss_ub_ne_n] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:directed_potential_loss_ub_ne_p] = Dict{Int, JuMP.ConstraintRef}()
    end

    L = wm.ref[:nw][n][:links][a]["length"]
    resistances = wm.ref[:nw][n][:resistance][a]

    dhn = wm.var[:nw][n][:dhn][a]
    qn_ne_ub = JuMP.upper_bound.(wm.var[:nw][n][:qn_ne][a])
    slopes_n = resistances .* qn_ne_ub.^(alpha - 1.0)
    rhs_n = sum(slopes_n .* wm.var[:nw][n][:qn_ne][a])
    con_n = JuMP.@constraint(wm.model, inv(L) * dhn - rhs_n <= 0.0)
    wm.con[:nw][n][:directed_potential_loss_ub_ne_n] = con_n

    dhp = wm.var[:nw][n][:dhp][a]
    qp_ne_ub = JuMP.upper_bound.(wm.var[:nw][n][:qp_ne][a])
    slopes_p = resistances .* qp_ne_ub.^(alpha - 1.0)
    rhs_p = sum(slopes_p .* wm.var[:nw][n][:qp_ne][a])
    con_p = JuMP.@constraint(wm.model, inv(L) * dhp - rhs_p <= 0.0)
    wm.con[:nw][n][:directed_potential_loss_ub_ne_p] = con_p
end

function constraint_directed_potential_loss_ub(wm::GenericWaterModel, a::Int, n::Int=wm.cnw; alpha::Float64=1.852)
    if !haskey(wm.con[:nw][n], :directed_potential_loss_ub_n)
        wm.con[:nw][n][:directed_potential_loss_ub_n] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:directed_potential_loss_ub_p] = Dict{Int, JuMP.ConstraintRef}()
    end

    L = wm.ref[:nw][n][:links][a]["length"]
    r = maximum(wm.ref[:nw][n][:resistance][a])

    dhn = wm.var[:nw][n][:dhn][a]
    qn_ub = JuMP.upper_bound(wm.var[:nw][n][:qn][a])
    rhs_n = r * qn_ub^(alpha - 1.0) * wm.var[:nw][n][:qn][a]
    con_n = JuMP.@constraint(wm.model, inv(L) * dhn - rhs_n <= 0.0)
    wm.con[:nw][n][:directed_potential_loss_ub_n] = con_n

    dhp = wm.var[:nw][n][:dhp][a]
    qp_ub = JuMP.upper_bound(wm.var[:nw][n][:qp][a])
    rhs_p = r * qp_ub^(alpha - 1.0) * wm.var[:nw][n][:qp][a]
    con_p = JuMP.@constraint(wm.model, inv(L) * dhp - rhs_p <= 0.0)
    wm.con[:nw][n][:directed_potential_loss_ub_p] = con_p
end

function constraint_link_undirected_flow_ne(wm::GenericWaterModel, a::Int, n::Int=wm.cnw; alpha::Float64=1.852)
    if !haskey(wm.con[:nw][n], :link_undirected_flow_ne)
        wm.con[:nw][n][:link_undirected_flow_ne] = Dict{Int, JuMP.ConstraintRef}()
    end

    q_ne = wm.var[:nw][n][:q_ne][a]
    q = wm.var[:nw][n][:q][a]
    con = JuMP.@constraint(wm.model, sum(q_ne) - q == 0.0)
    wm.con[:nw][n][:link_undirected_flow_ne] = con
end

function constraint_link_directed_flow_ne(wm::GenericWaterModel, a::Int, n::Int=wm.cnw; alpha::Float64=1.852)
    if !haskey(wm.con[:nw][n], :link_directed_flow_n_ne)
        wm.con[:nw][n][:link_directed_flow_n_ne] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:link_directed_flow_p_ne] = Dict{Int, JuMP.ConstraintRef}()
    end

    qn_ne = wm.var[:nw][n][:qn_ne][a]
    qn = wm.var[:nw][n][:qn][a]
    con_n = JuMP.@constraint(wm.model, sum(qn_ne) - qn == 0.0)
    wm.con[:nw][n][:link_directed_flow_n_ne] = con_n

    qp_ne = wm.var[:nw][n][:qp_ne][a]
    qp = wm.var[:nw][n][:qp][a]
    con_p = JuMP.@constraint(wm.model, sum(qp_ne) - qp == 0.0)
    wm.con[:nw][n][:link_directed_flow_p_ne] = con_p
end

function constraint_link_directed_flow(wm::GenericWaterModel, a::Int, n::Int=wm.cnw; alpha::Float64=1.852)
    if !haskey(wm.con[:nw][n], :link_directed_flow)
        wm.con[:nw][n][:link_directed_flow] = Dict{Int, JuMP.ConstraintRef}()
    end

    q = wm.var[:nw][n][:q][a]
    qn = wm.var[:nw][n][:qn][a]
    qp = wm.var[:nw][n][:qp][a]
    con = JuMP.@constraint(wm.model, (qp - qn) - q == 0.0)
    wm.con[:nw][n][:link_directed_flow] = con
end

function constraint_undirected_sink_flow(wm::GenericWaterModel, i::Int, n::Int=wm.cnw)
end

function constraint_undirected_source_flow(wm::GenericWaterModel, i::Int, n::Int=wm.cnw)
end

"Constraint to ensure at least one direction is set to take flow away from a source."
function constraint_directed_source_flow(wm::GenericWaterModel, i::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :directed_source_flow)
        wm.con[:nw][n][:directed_source_flow] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Collect the required variables.
    x_dir = wm.var[:nw][n][:x_dir]
    out_arcs = filter(a -> i == a.second["f_id"], wm.ref[:nw][n][:links])
    out = Array{JuMP.VariableRef}([x_dir[a] for a in keys(out_arcs)])
    in_arcs = filter(a -> i == a.second["t_id"], wm.ref[:nw][n][:links])
    in = Array{JuMP.VariableRef}([x_dir[a] for a in keys(in_arcs)])

    # Add the source flow direction constraint.
    con = JuMP.@constraint(wm.model, sum(out) - sum(in) >= 1.0 - length(in))
    wm.con[:nw][n][:directed_source_flow][i] = con
end

"Constraint to ensure at least one direction is set to take flow to a junction with demand."
function constraint_directed_sink_flow(wm::GenericWaterModel, i::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :directed_sink_flow)
        wm.con[:nw][n][:directed_sink_flow] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Collect the required variables.
    x_dir = wm.var[:nw][n][:x_dir]
    out_arcs = filter(a -> i == a.second["f_id"], wm.ref[:nw][n][:links])
    out = Array{JuMP.VariableRef}([x_dir[a] for a in keys(out_arcs)])
    in_arcs = filter(a -> i == a.second["t_id"], wm.ref[:nw][n][:links])
    in = Array{JuMP.VariableRef}([x_dir[a] for a in keys(in_arcs)])

    # Add the sink flow direction constraint.
    con = JuMP.@constraint(wm.model, sum(in) - sum(out) >= 1.0 - length(out))
    wm.con[:nw][n][:directed_sink_flow][i] = con
end
