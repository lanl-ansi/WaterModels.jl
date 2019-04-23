########################################################################
# This file defines commonly-used constraints for water systems models.
########################################################################

function constraint_directed_flow_conservation(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Initialize the flow sum expression.
    flow_sum = JuMP.AffExpr(0.0)

    # Add all incoming flow to node i.
    for (a, conn) in filter(is_in_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, -wm.var[:nw][n][:q⁻][a])
        JuMP.add_to_expression!(flow_sum, wm.var[:nw][n][:q⁺][a])
    end

    # Subtract all outgoing flow from node i.
    for (a, conn) in filter(is_out_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, wm.var[:nw][n][:q⁻][a])
        JuMP.add_to_expression!(flow_sum, -wm.var[:nw][n][:q⁺][a])
    end

    # Create the constraint dictionary if necessary.
    if !haskey(wm.con[:nw][n], :flow_conservation)
        wm.con[:nw][n][:flow_conservation] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Add the flow conservation constraint.
    demand = wm.ref[:nw][n][:junctions][i]["demand"]
    con = JuMP.@constraint(wm.model, flow_sum == demand)
    wm.con[:nw][n][:flow_conservation][i] = con
end

function constraint_directed_flow_conservation_ne(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Initialize the flow sum expression.
    flow_sum = JuMP.AffExpr(0.0)

    # Add all incoming flow to node i.
    for (a, conn) in filter(is_in_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, -sum(wm.var[:nw][n][:q_ne⁻][a, :]))
        JuMP.add_to_expression!(flow_sum, sum(wm.var[:nw][n][:q_ne⁺][a, :]))
    end

    # Subtract all outgoing flow from node i.
    for (a, conn) in filter(is_out_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, sum(wm.var[:nw][n][:q_ne⁻][a, :]))
        JuMP.add_to_expression!(flow_sum, -sum(wm.var[:nw][n][:q_ne⁺][a, :]))
    end

    # Create the constraint dictionary if necessary.
    if !haskey(wm.con[:nw][n], :flow_conservation)
        wm.con[:nw][n][:flow_conservation] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Add the flow conservation constraint.
    demand = wm.ref[:nw][n][:junctions][i]["demand"]
    con = JuMP.@constraint(wm.model, flow_sum == demand)
    wm.con[:nw][n][:flow_conservation][i] = con
end

function constraint_resistance_selection(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    if !haskey(wm.con[:nw][n], :resistance_selection_sum)
        wm.con[:nw][n][:resistance_selection_sum] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:resistance_selection_n] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        wm.con[:nw][n][:resistance_selection_p] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    wm.con[:nw][n][:resistance_selection_n][a] = Dict{Int, JuMP.ConstraintRef}()
    wm.con[:nw][n][:resistance_selection_p][a] = Dict{Int, JuMP.ConstraintRef}()

    con_sum = JuMP.@constraint(wm.model, sum(wm.var[:nw][n][:xr][a]) == 1)
    wm.con[:nw][n][:resistance_selection_sum][a] = con_sum

    for r in 1:length(wm.ref[:nw][n][:resistance][a])
        x_a_res = wm.var[:nw][n][:xr][a][r]

        q_n_a_r = wm.var[:nw][n][:q_ne⁻][a, r]
        q_n_a_r_ub = JuMP.upper_bound(q_n_a_r)
        con_n = JuMP.@constraint(wm.model, q_n_a_r <= q_n_a_r_ub * x_a_res)
        wm.con[:nw][n][:resistance_selection_n][a][r] = con_n

        q_p_a_r = wm.var[:nw][n][:q_ne⁺][a, r]
        q_p_a_r_ub = JuMP.upper_bound(q_p_a_r)
        con_p = JuMP.@constraint(wm.model, q_p_a_r <= q_p_a_r_ub * x_a_res)
        wm.con[:nw][n][:resistance_selection_p][a][r] = con_p
    end
end

function constraint_flow_direction_selection(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :flow_direction_selection_n)
        wm.con[:nw][n][:flow_direction_selection_n] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        wm.con[:nw][n][:flow_direction_selection_p] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    wm.con[:nw][n][:flow_direction_selection_n][a] = Dict{Int, JuMP.ConstraintRef}()
    wm.con[:nw][n][:flow_direction_selection_p][a] = Dict{Int, JuMP.ConstraintRef}()

    for r in 1:length(wm.ref[:nw][n][:resistance][a])
        x_a_dir = wm.var[:nw][n][:dir][a]

        q_n_a_r = wm.var[:nw][n][:q⁻][a][r]
        q_n_a_r_ub = JuMP.upper_bound(q_n_a_r)
        con_n = JuMP.@constraint(wm.model, q_n_a_r <= q_n_a_r_ub * (1 - x_a_dir))
        wm.con[:nw][n][:flow_direction_selection_n][a][r] = con_n

        q_p_a_r = wm.var[:nw][n][:q⁺][a][r]
        q_p_a_r_ub = JuMP.upper_bound(q_p_a_r)
        con_p = JuMP.@constraint(wm.model, q_p_a_r <= q_p_a_r_ub * x_a_dir)
        wm.con[:nw][n][:flow_direction_selection_p][a][r] = con_p
    end
end
