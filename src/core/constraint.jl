########################################################################
# This file defines commonly-used constraints for water systems models.
########################################################################

function constraint_undirected_flow_conservation(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Create the constraint dictionary if necessary.
    if !haskey(wm.con[:nw][n], :flow_conservation)
        con(wm, n)[:flow_conservation] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Initialize the flow sum expression.
    flow_sum = JuMP.AffExpr(0.0)

    # Add all incoming flow to node i.
    for (a, link) in filter(is_in_node(i), ref(wm, n, :links))
        JuMP.add_to_expression!(flow_sum, var(wm, n, :q, a))
    end

    # Subtract all outgoing flow from node i.
    for (a, link) in filter(is_out_node(i), ref(wm, n, :links))
        JuMP.add_to_expression!(flow_sum, -var(wm, n, :q, a))
    end

    # Add the flow conservation constraint.
    demand = ref(wm, n, :junctions, i)["base_demand"]
    c = JuMP.@constraint(wm.model, flow_sum == demand)
    con(wm, n, :flow_conservation)[i] = c
end

function constraint_undirected_flow_conservation_ne(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Create the constraint dictionary if necessary.
    if !haskey(wm.con[:nw][n], :flow_conservation_ne)
        con(wm, n)[:flow_conservation_ne] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Initialize the flow sum expression.
    flow_sum = JuMP.AffExpr(0.0)

    # Add all incoming flow to node i.
    for (a, conn) in filter(is_in_node(i), ref(wm, n, :links))
        JuMP.add_to_expression!(flow_sum, sum(var(wm, n, :q_ne, a)))
    end

    # Subtract all outgoing flow from node i.
    for (a, conn) in filter(is_out_node(i), ref(wm, n, :links))
        JuMP.add_to_expression!(flow_sum, -sum(var(wm, n, :q_ne, a)))
    end

    # Add the flow conservation constraint.
    demand = ref(wm, n, :junctions, i)["base_demand"]
    c = JuMP.@constraint(wm.model, flow_sum == demand)
    con(wm, n, :flow_conservation_ne)[i] = c
end

function constraint_directed_flow_conservation(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Create the constraint dictionary if necessary.
    if !haskey(wm.con[:nw][n], :flow_conservation)
        con(wm, n)[:flow_conservation] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Initialize the flow sum expression.
    flow_sum = JuMP.AffExpr(0.0)

    # Add all incoming flow to node i.
    for (a, conn) in filter(is_in_node(i), ref(wm, n, :links))
        JuMP.add_to_expression!(flow_sum, var(wm, n, :qp, a))
        JuMP.add_to_expression!(flow_sum, -var(wm, n, :qn, a))
    end

    # Subtract all outgoing flow from node i.
    for (a, conn) in filter(is_out_node(i), ref(wm, n, :links))
        JuMP.add_to_expression!(flow_sum, -var(wm, n, :qp, a))
        JuMP.add_to_expression!(flow_sum, var(wm, n, :qn, a))
    end

    # Add the flow conservation constraint.
    demand = ref(wm, n, :junctions, i)["base_demand"]
    c = JuMP.@constraint(wm.model, flow_sum == demand)
    con(wm, n, :flow_conservation)[i] = c
end

function constraint_directed_flow_conservation_ne(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Create the constraint dictionary if necessary.
    if !haskey(wm.con[:nw][n], :flow_conservation_ne)
        con(wm, n)[:flow_conservation_ne] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Initialize the flow sum expression.
    flow_sum = JuMP.AffExpr(0.0)

    # Add all incoming flow to node i.
    for (a, conn) in filter(is_in_node(i), ref(wm, n, :links))
        JuMP.add_to_expression!(flow_sum, sum(var(wm, n, :qp_ne, a)))
        JuMP.add_to_expression!(flow_sum, -sum(var(wm, n, :qn_ne, a)))
    end

    # Subtract all outgoing flow from node i.
    for (a, conn) in filter(is_out_node(i), ref(wm, n, :links))
        JuMP.add_to_expression!(flow_sum, -sum(var(wm, n, :qp_ne, a)))
        JuMP.add_to_expression!(flow_sum, sum(var(wm, n, :qn_ne, a)))
    end

    # Add the flow conservation constraint.
    demand = ref(wm, n, :junctions, i)["base_demand"]
    c = JuMP.@constraint(wm.model, flow_sum == demand)
    con(wm, n, :flow_conservation_ne)[i] = c
end

function constraint_directed_resistance_selection_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    if !haskey(wm.con[:nw][n], :resistance_selection_sum)
        con(wm, n)[:resistance_selection_sum] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:resistance_selection_p] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        con(wm, n)[:resistance_selection_n] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    con(wm, n, :resistance_selection_p)[a] = Dict{Int, JuMP.ConstraintRef}()
    con(wm, n, :resistance_selection_n)[a] = Dict{Int, JuMP.ConstraintRef}()

    con_sum = JuMP.@constraint(wm.model, sum(var(wm, n, :x_res, a)) == 1.0)
    con(wm, n, :resistance_selection_sum)[a] = con_sum

    for r in 1:length(ref(wm, n, :resistance, a))
        x_res = var(wm, n, :x_res, a)[r]

        qp_ne = var(wm, n, :qp_ne, a)[r]
        qp_ne_ub = JuMP.upper_bound(qp_ne)
        con_p = JuMP.@constraint(wm.model, qp_ne - qp_ne_ub * x_res <= 0.0)
        con(wm, n, :resistance_selection_p, a)[r] = con_p

        qn_ne = var(wm, n, :qn_ne, a)[r]
        qn_ne_ub = JuMP.upper_bound(qn_ne)
        con_n = JuMP.@constraint(wm.model, qn_ne - qn_ne_ub * x_res <= 0.0)
        con(wm, n, :resistance_selection_n, a)[r] = con_n
    end
end

function constraint_undirected_resistance_selection_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    if !haskey(wm.con[:nw][n], :resistance_selection_sum)
        con(wm, n)[:resistance_selection_sum] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:resistance_selection_lb] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        con(wm, n)[:resistance_selection_ub] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    con_sum = JuMP.@constraint(wm.model, sum(var(wm, n, :x_res, a)) == 1.0)
    con(wm, n, :resistance_selection_sum)[a] = con_sum

    con(wm, n, :resistance_selection_lb)[a] = Dict{Int, JuMP.ConstraintRef}()
    con(wm, n, :resistance_selection_ub)[a] = Dict{Int, JuMP.ConstraintRef}()

    for r in 1:length(ref(wm, n, :resistance, a))
        x_res = var(wm, n, :x_res, a)[r]

        q_ne = var(wm, n, :q_ne, a)[r]
        q_ne_lb = JuMP.lower_bound(q_ne)
        con_lb = JuMP.@constraint(wm.model, q_ne - q_ne_lb * x_res >= 0.0)
        con(wm, n, :resistance_selection_lb, a)[r] = con_lb

        q_ne = var(wm, n, :q_ne, a)[r]
        q_ne_ub = JuMP.upper_bound(q_ne)
        con_ub = JuMP.@constraint(wm.model, q_ne - q_ne_ub * x_res <= 0.0)
        con(wm, n, :resistance_selection_ub, a)[r] = con_ub
    end
end

function constraint_flow_direction_selection(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :flow_direction_selection_n)
        con(wm, n)[:flow_direction_selection_n] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:flow_direction_selection_p] = Dict{Int, JuMP.ConstraintRef}()
    end

    x_dir = var(wm, n, :x_dir, a)

    qp = var(wm, n, :qp, a)
    qp_ub = JuMP.upper_bound(qp)
    con_p = JuMP.@constraint(wm.model, qp - qp_ub * x_dir <= 0.0)
    con(wm, n, :flow_direction_selection_p)[a] = con_p

    qn = var(wm, n, :qn, a)
    qn_ub = JuMP.upper_bound(qn)
    con_n = JuMP.@constraint(wm.model, qn - qn_ub * (1.0 - x_dir) <= 0.0)
    con(wm, n, :flow_direction_selection_n)[a] = con_n
end

function constraint_flow_direction_selection_ne(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :flow_direction_selection_ne_n)
        con(wm, n)[:flow_direction_selection_ne_p] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        con(wm, n)[:flow_direction_selection_ne_n] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    con(wm, n, :flow_direction_selection_ne_p)[a] = Dict{Int, JuMP.ConstraintRef}()
    con(wm, n, :flow_direction_selection_ne_n)[a] = Dict{Int, JuMP.ConstraintRef}()

    for r in 1:length(ref(wm, n, :resistance, a))
        x_dir = var(wm, n, :x_dir, a)

        qp_ne = var(wm, n, :qp_ne, a)[r]
        qp_ne_ub = JuMP.upper_bound(qp_ne)
        con_p = JuMP.@constraint(wm.model, qp_ne - qp_ne_ub * x_dir <= 0.0)
        con(wm, n, :flow_direction_selection_ne_p, a)[r] = con_p

        qn_ne = var(wm, n, :qn_ne, a)[r]
        qn_ne_ub = JuMP.upper_bound(qn_ne)
        con_n = JuMP.@constraint(wm.model, qn_ne - qn_ne_ub * (1.0 - x_dir) <= 0.0)
        con(wm, n, :flow_direction_selection_ne_n, a)[r] = con_n
    end
end

function constraint_directed_head_difference(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :head_difference_1)
        con(wm, n)[:head_difference_1] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:head_difference_2] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:head_difference_3] = Dict{Int, JuMP.ConstraintRef}()
    end

    i = ref(wm, n, :links, a)["f_id"]

    if i in collect(ids(wm, n, :reservoirs))
        h_i = ref(wm, n, :reservoirs, i)["base_head"]
    else
        h_i = var(wm, n, :h, i)
    end

    j = ref(wm, n, :links, a)["t_id"]

    if j in collect(ids(wm, n, :reservoirs))
        h_j = ref(wm, n, :reservoirs, j)["base_head"]
    else
        h_j = var(wm, n, :h, j)
    end

    x_dir = var(wm, n, :x_dir, a)

    dhp = var(wm, n, :dhp, a)
    dhp_ub = JuMP.upper_bound(dhp)
    con_1 = JuMP.@constraint(wm.model, dhp - dhp_ub * x_dir <= 0.0)
    con(wm, n, :head_difference_1)[a] = con_1

    dhn = var(wm, n, :dhn, a)
    dhn_ub = JuMP.upper_bound(dhn)
    con_2 = JuMP.@constraint(wm.model, dhn - dhn_ub * (1.0 - x_dir) <= 0.0)
    con(wm, n, :head_difference_2)[a] = con_2

    con_3 = JuMP.@constraint(wm.model, (h_i - h_j) - (dhp - dhn) == 0.0)
    con(wm, n, :head_difference_3)[a] = con_3
end

function constraint_directed_potential_loss_ub_ne(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :directed_potential_loss_ub_ne_n)
        con(wm, n)[:directed_potential_loss_ub_ne_n] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:directed_potential_loss_ub_ne_p] = Dict{Int, JuMP.ConstraintRef}()
    end

    alpha = ref(wm, n, :alpha)
    L = ref(wm, n, :links, a)["length"]
    resistances = ref(wm, n, :resistance, a)

    dhp = var(wm, n, :dhp, a)
    qp_ne_ub = JuMP.upper_bound.(var(wm, n, :qp_ne, a))
    slopes_p = resistances .* qp_ne_ub.^(alpha - 1.0)
    rhs_p = sum(slopes_p .* var(wm, n, :qp_ne, a))
    con_p = JuMP.@constraint(wm.model, inv(L) * dhp - rhs_p <= 0.0)
    con(wm, n, :directed_potential_loss_ub_ne_p)[a] = con_p

    dhn = var(wm, n, :dhn, a)
    qn_ne_ub = JuMP.upper_bound.(var(wm, n, :qn_ne, a))
    slopes_n = resistances .* qn_ne_ub.^(alpha - 1.0)
    rhs_n = sum(slopes_n .* var(wm, n, :qn_ne, a))
    con_n = JuMP.@constraint(wm.model, inv(L) * dhn - rhs_n <= 0.0)
    con(wm, n, :directed_potential_loss_ub_ne_n)[a] = con_n
end

function constraint_directed_potential_loss_ub(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :directed_potential_loss_ub_n)
        con(wm, n)[:directed_potential_loss_ub_p] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:directed_potential_loss_ub_n] = Dict{Int, JuMP.ConstraintRef}()
    end

    alpha = ref(wm, n, :alpha)
    L = ref(wm, n, :links, a)["length"]
    r = maximum(ref(wm, n, :resistance, a))

    dhp = var(wm, n, :dhp, a)
    qp_ub = JuMP.upper_bound(var(wm, n, :qp, a))
    rhs_p = r * qp_ub^(alpha - 1.0) * var(wm, n, :qp, a)
    con_p = JuMP.@constraint(wm.model, inv(L) * dhp - rhs_p <= 0.0)
    con(wm, n, :directed_potential_loss_ub_p)[a] = con_p

    dhn = var(wm, n, :dhn, a)
    qn_ub = JuMP.upper_bound(var(wm, n, :qn, a))
    rhs_n = r * qn_ub^(alpha - 1.0) * var(wm, n, :qn, a)
    con_n = JuMP.@constraint(wm.model, inv(L) * dhn - rhs_n <= 0.0)
    con(wm, n, :directed_potential_loss_ub_n)[a] = con_n
end

function constraint_link_undirected_flow_ne(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :link_undirected_flow_ne)
        con(wm, n)[:link_undirected_flow_ne] = Dict{Int, JuMP.ConstraintRef}()
    end

    q_ne = var(wm, n, :q_ne, a)
    q = var(wm, n, :q, a)
    c = JuMP.@constraint(wm.model, sum(q_ne) - q == 0.0)
    con(wm, n, :link_undirected_flow_ne)[a] = c
end

function constraint_link_directed_flow_ne(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :link_directed_flow_n_ne)
        con(wm, n)[:link_directed_flow_p_ne] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:link_directed_flow_n_ne] = Dict{Int, JuMP.ConstraintRef}()
    end

    qp_ne = var(wm, n, :qp_ne, a)
    qp = var(wm, n, :qp, a)
    con_p = JuMP.@constraint(wm.model, sum(qp_ne) - qp == 0.0)
    con(wm, n, :link_directed_flow_p_ne)[a] = con_p

    qn_ne = var(wm, n, :qn_ne, a)
    qn = var(wm, n, :qn, a)
    con_n = JuMP.@constraint(wm.model, sum(qn_ne) - qn == 0.0)
    con(wm, n, :link_directed_flow_n_ne)[a] = con_n
end

function constraint_link_directed_flow(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :link_directed_flow)
        con(wm, n)[:link_directed_flow] = Dict{Int, JuMP.ConstraintRef}()
    end

    q = var(wm, n, :q, a)
    qp = var(wm, n, :qp, a)
    qn = var(wm, n, :qn, a)
    c = JuMP.@constraint(wm.model, (qp - qn) - q == 0.0)
    con(wm, n, :link_directed_flow)[a] = c
end

function constraint_undirected_sink_flow(wm::GenericWaterModel, i::Int, n::Int=wm.cnw)
end

function constraint_undirected_source_flow(wm::GenericWaterModel, i::Int, n::Int=wm.cnw)
end

"Constraint to ensure at least one direction is set to take flow away from a source."
function constraint_directed_source_flow(wm::GenericWaterModel, i::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :directed_source_flow)
        con(wm, n)[:directed_source_flow] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Collect the required variables.
    x_dir = var(wm, n, :x_dir)
    out_arcs = filter(a -> i == a.second["f_id"], ref(wm, n, :links))
    out = Array{JuMP.VariableRef}([x_dir[a] for a in keys(out_arcs)])
    in_arcs = filter(a -> i == a.second["t_id"], ref(wm, n, :links))
    in = Array{JuMP.VariableRef}([x_dir[a] for a in keys(in_arcs)])

    # Add the source flow direction constraint.
    c = JuMP.@constraint(wm.model, sum(out) - sum(in) >= 1.0 - length(in))
    con(wm, n, :directed_source_flow)[i] = c
end

"Constraint to ensure at least one direction is set to take flow to a junction with demand."
function constraint_directed_sink_flow(wm::GenericWaterModel, i::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :directed_sink_flow)
        con(wm, n)[:directed_sink_flow] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Collect the required variables.
    x_dir = var(wm, n, :x_dir)
    out_arcs = filter(a -> i == a.second["f_id"], ref(wm, n, :links))
    out = Array{JuMP.VariableRef}([x_dir[a] for a in keys(out_arcs)])
    in_arcs = filter(a -> i == a.second["t_id"], ref(wm, n, :links))
    in = Array{JuMP.VariableRef}([x_dir[a] for a in keys(in_arcs)])

    # Add the sink flow direction constraint.
    c = JuMP.@constraint(wm.model, sum(in) - sum(out) >= 1.0 - length(out))
    con(wm, n, :directed_sink_flow)[i] = c
end
