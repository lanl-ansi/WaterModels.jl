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
    for (a, conn) in filter(is_in_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, wm.var[:nw][n][:q][a])
    end

    # Subtract all outgoing flow from node i.
    for (a, conn) in filter(is_out_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, -wm.var[:nw][n][:q][a])
    end

    # Add the flow conservation constraint.
    demand = wm.ref[:nw][n][:junctions][i]["demand"]
    con = JuMP.@constraint(wm.model, flow_sum == demand)
    wm.con[:nw][n][:flow_conservation][i] = con
end

function constraint_undirected_flow_conservation_ne(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Create the constraint dictionary if necessary.
    if !haskey(wm.con[:nw][n], :flow_conservationⁿᵉ)
        wm.con[:nw][n][:flow_conservationⁿᵉ] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Initialize the flow sum expression.
    flow_sum = JuMP.AffExpr(0.0)

    # Add all incoming flow to node i.
    for (a, conn) in filter(is_in_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, sum(wm.var[:nw][n][:qⁿᵉ][a]))
    end

    # Subtract all outgoing flow from node i.
    for (a, conn) in filter(is_out_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, -sum(wm.var[:nw][n][:qⁿᵉ][a]))
    end

    # Add the flow conservation constraint.
    demand = wm.ref[:nw][n][:junctions][i]["demand"]
    con = JuMP.@constraint(wm.model, flow_sum == demand)
    wm.con[:nw][n][:flow_conservationⁿᵉ][i] = con
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
        JuMP.add_to_expression!(flow_sum, -wm.var[:nw][n][:q⁻][a])
        JuMP.add_to_expression!(flow_sum, wm.var[:nw][n][:q⁺][a])
    end

    # Subtract all outgoing flow from node i.
    for (a, conn) in filter(is_out_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, wm.var[:nw][n][:q⁻][a])
        JuMP.add_to_expression!(flow_sum, -wm.var[:nw][n][:q⁺][a])
    end

    # Add the flow conservation constraint.
    demand = wm.ref[:nw][n][:junctions][i]["demand"]
    con = JuMP.@constraint(wm.model, flow_sum == demand)
    wm.con[:nw][n][:flow_conservation][i] = con
end

function constraint_directed_flow_conservation_ne(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Create the constraint dictionary if necessary.
    if !haskey(wm.con[:nw][n], :flow_conservationⁿᵉ)
        wm.con[:nw][n][:flow_conservationⁿᵉ] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Initialize the flow sum expression.
    flow_sum = JuMP.AffExpr(0.0)

    # Add all incoming flow to node i.
    for (a, conn) in filter(is_in_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, -sum(wm.var[:nw][n][:qⁿᵉ⁻][a]))
        JuMP.add_to_expression!(flow_sum, sum(wm.var[:nw][n][:qⁿᵉ⁺][a]))
    end

    # Subtract all outgoing flow from node i.
    for (a, conn) in filter(is_out_node(i), wm.ref[:nw][n][:links])
        JuMP.add_to_expression!(flow_sum, sum(wm.var[:nw][n][:qⁿᵉ⁻][a]))
        JuMP.add_to_expression!(flow_sum, -sum(wm.var[:nw][n][:qⁿᵉ⁺][a]))
    end

    # Add the flow conservation constraint.
    demand = wm.ref[:nw][n][:junctions][i]["demand"]
    con = JuMP.@constraint(wm.model, flow_sum == demand)
    wm.con[:nw][n][:flow_conservationⁿᵉ][i] = con
end

function constraint_directed_resistance_selection_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    if !haskey(wm.con[:nw][n], :resistance_selection_sum)
        wm.con[:nw][n][:resistance_selection_sum] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:resistance_selection⁻] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        wm.con[:nw][n][:resistance_selection⁺] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    wm.con[:nw][n][:resistance_selection⁻][a] = Dict{Int, JuMP.ConstraintRef}()
    wm.con[:nw][n][:resistance_selection⁺][a] = Dict{Int, JuMP.ConstraintRef}()

    con_sum = JuMP.@constraint(wm.model, sum(wm.var[:nw][n][:xʳᵉˢ][a]) == 1.0)
    wm.con[:nw][n][:resistance_selection_sum][a] = con_sum

    for r in 1:length(wm.ref[:nw][n][:resistance][a])
        xʳᵉˢ = wm.var[:nw][n][:xʳᵉˢ][a][r]

        qⁿᵉ⁻ = wm.var[:nw][n][:qⁿᵉ⁻][a][r]
        q̅ⁿᵉ⁻ = JuMP.upper_bound(qⁿᵉ⁻)
        con⁻ = JuMP.@constraint(wm.model, qⁿᵉ⁻ - q̅ⁿᵉ⁻ * xʳᵉˢ <= 0.0)
        wm.con[:nw][n][:resistance_selection⁻][a][r] = con⁻

        qⁿᵉ⁺ = wm.var[:nw][n][:qⁿᵉ⁺][a][r]
        q̅ⁿᵉ⁺ = JuMP.upper_bound(qⁿᵉ⁺)
        con⁺ = JuMP.@constraint(wm.model, qⁿᵉ⁺ - q̅ⁿᵉ⁺ * xʳᵉˢ <= 0.0)
        wm.con[:nw][n][:resistance_selection⁺][a][r] = con⁺
    end
end

function constraint_undirected_resistance_selection_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    if !haskey(wm.con[:nw][n], :resistance_selection_sum)
        wm.con[:nw][n][:resistance_selection_sum] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:resistance_selection_lb] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        wm.con[:nw][n][:resistance_selection_ub] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    con_sum = JuMP.@constraint(wm.model, sum(wm.var[:nw][n][:xʳᵉˢ][a]) == 1.0)
    wm.con[:nw][n][:resistance_selection_sum][a] = con_sum

    wm.con[:nw][n][:resistance_selection_lb][a] = Dict{Int, JuMP.ConstraintRef}()
    wm.con[:nw][n][:resistance_selection_ub][a] = Dict{Int, JuMP.ConstraintRef}()

    for r in 1:length(wm.ref[:nw][n][:resistance][a])
        xʳᵉˢ = wm.var[:nw][n][:xʳᵉˢ][a][r]

        qⁿᵉ = wm.var[:nw][n][:qⁿᵉ][a][r]
        q͟ⁿᵉ = JuMP.lower_bound(qⁿᵉ)
        con_lb = JuMP.@constraint(wm.model, qⁿᵉ - q͟ⁿᵉ * xʳᵉˢ >= 0.0)
        wm.con[:nw][n][:resistance_selection_lb][a][r] = con_lb

        qⁿᵉ = wm.var[:nw][n][:qⁿᵉ][a][r]
        q̅ⁿᵉ = JuMP.upper_bound(qⁿᵉ)
        con_ub = JuMP.@constraint(wm.model, qⁿᵉ - q̅ⁿᵉ * xʳᵉˢ <= 0.0)
        wm.con[:nw][n][:resistance_selection_ub][a][r] = con_ub
    end
end

function constraint_flow_direction_selection(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :flow_direction_selection⁻)
        wm.con[:nw][n][:flow_direction_selection⁻] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:flow_direction_selection⁺] = Dict{Int, JuMP.ConstraintRef}()
    end

    xᵈᶦʳ = wm.var[:nw][n][:xᵈᶦʳ][a]

    q⁻ = wm.var[:nw][n][:q⁻][a]
    q̅⁻ = JuMP.upper_bound(q⁻)
    con⁻ = JuMP.@constraint(wm.model, q⁻ - q̅⁻ * (1.0 - xᵈᶦʳ) <= 0.0)
    wm.con[:nw][n][:flow_direction_selection⁻][a] = con⁻

    q⁺ = wm.var[:nw][n][:q⁺][a]
    q̅⁺ = JuMP.upper_bound(q⁺)
    con⁺ = JuMP.@constraint(wm.model, q⁺ - q̅⁺ * xᵈᶦʳ <= 0.0)
    wm.con[:nw][n][:flow_direction_selection⁺][a] = con⁺
end

function constraint_flow_direction_selection_ne(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :flow_direction_selection_ne⁻)
        wm.con[:nw][n][:flow_direction_selection_ne⁻] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        wm.con[:nw][n][:flow_direction_selection_ne⁺] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    wm.con[:nw][n][:flow_direction_selection_ne⁻][a] = Dict{Int, JuMP.ConstraintRef}()
    wm.con[:nw][n][:flow_direction_selection_ne⁺][a] = Dict{Int, JuMP.ConstraintRef}()

    for r in 1:length(wm.ref[:nw][n][:resistance][a])
        xᵈᶦʳ = wm.var[:nw][n][:xᵈᶦʳ][a]

        qⁿᵉ⁻ = wm.var[:nw][n][:qⁿᵉ⁻][a][r]
        q̅ⁿᵉ⁻ = JuMP.upper_bound(qⁿᵉ⁻)
        con⁻ = JuMP.@constraint(wm.model, qⁿᵉ⁻ - q̅ⁿᵉ⁻ * (1.0 - xᵈᶦʳ) <= 0.0)
        wm.con[:nw][n][:flow_direction_selection_ne⁻][a][r] = con⁻

        qⁿᵉ⁺ = wm.var[:nw][n][:qⁿᵉ⁺][a][r]
        q̅ⁿᵉ⁺ = JuMP.upper_bound(qⁿᵉ⁺)
        con⁺ = JuMP.@constraint(wm.model, qⁿᵉ⁺ - q̅ⁿᵉ⁺ * xᵈᶦʳ <= 0.0)
        wm.con[:nw][n][:flow_direction_selection_ne⁺][a][r] = con⁺
    end
end

function constraint_directed_head_difference(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(wm.con[:nw][n], :head_difference_1)
        wm.con[:nw][n][:head_difference_1] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:head_difference_2] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:head_difference_3] = Dict{Int, JuMP.ConstraintRef}()
    end

    i = wm.ref[:nw][n][:links][a]["node1"]

    if i in collect(ids(wm, n, :reservoirs))
        hᵢ = wm.ref[:nw][n][:reservoirs][i]["head"]
    else
        hᵢ = wm.var[:nw][n][:h][i]
    end

    j = wm.ref[:nw][n][:links][a]["node2"]

    if j in collect(ids(wm, n, :reservoirs))
        hⱼ = wm.ref[:nw][n][:reservoirs][j]["head"]
    else
        hⱼ = wm.var[:nw][n][:h][j]
    end

    xᵈᶦʳ = wm.var[:nw][n][:xᵈᶦʳ][a]

    Δh⁻ = wm.var[:nw][n][:Δh⁻][a]
    Δh̅⁻ = JuMP.upper_bound(Δh⁻)
    con_2 = JuMP.@constraint(wm.model, Δh⁻ - Δh̅⁻ * (1.0 - xᵈᶦʳ) <= 0.0)
    wm.con[:nw][n][:head_difference_2][a] = con_2

    Δh⁺ = wm.var[:nw][n][:Δh⁺][a]
    Δh̅⁺ = JuMP.upper_bound(Δh⁺)
    con_1 = JuMP.@constraint(wm.model, Δh⁺ - Δh̅⁺ * xᵈᶦʳ <= 0.0)
    wm.con[:nw][n][:head_difference_1][a] = con_1

    con_3 = JuMP.@constraint(wm.model, (hᵢ - hⱼ) - (Δh⁺ - Δh⁻) == 0.0)
    wm.con[:nw][n][:head_difference_3][a] = con_3
end

function constraint_directed_potential_loss_ub_ne(wm::GenericWaterModel, a::Int, n::Int=wm.cnw; alpha::Float64=1.852)
    if !haskey(wm.con[:nw][n], :directed_potential_loss_ub_ne⁻)
        wm.con[:nw][n][:directed_potential_loss_ub_ne⁻] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:directed_potential_loss_ub_ne⁺] = Dict{Int, JuMP.ConstraintRef}()
    end

    L = wm.ref[:nw][n][:links][a]["length"]
    resistances = wm.ref[:nw][n][:resistance][a]

    Δh⁻ = wm.var[:nw][n][:Δh⁻][a]
    q̅ⁿᵉ⁻ = JuMP.upper_bound.(wm.var[:nw][n][:qⁿᵉ⁻][a])
    slopes⁻ = resistances .* q̅ⁿᵉ⁻.^(alpha - 1.0)
    rhs⁻ = sum(slopes⁻ .* wm.var[:nw][n][:qⁿᵉ⁻][a])
    con⁻ = JuMP.@constraint(wm.model, inv(L) * Δh⁻ - rhs⁻ <= 0.0)
    wm.con[:nw][n][:directed_potential_loss_ub_ne⁻] = con⁻

    Δh⁺ = wm.var[:nw][n][:Δh⁺][a]
    q̅ⁿᵉ⁺ = JuMP.upper_bound.(wm.var[:nw][n][:qⁿᵉ⁺][a])
    slopes⁺ = resistances .* q̅ⁿᵉ⁺.^(alpha - 1.0)
    rhs⁺ = sum(slopes⁺ .* wm.var[:nw][n][:qⁿᵉ⁺][a])
    con⁺ = JuMP.@constraint(wm.model, inv(L) * Δh⁺ - rhs⁺ <= 0.0)
    wm.con[:nw][n][:directed_potential_loss_ub_ne⁺] = con⁺
end

function constraint_directed_potential_loss_ub(wm::GenericWaterModel, a::Int, n::Int=wm.cnw; alpha::Float64=1.852)
    if !haskey(wm.con[:nw][n], :directed_potential_loss_ub⁻)
        wm.con[:nw][n][:directed_potential_loss_ub⁻] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:directed_potential_loss_ub⁺] = Dict{Int, JuMP.ConstraintRef}()
    end

    L = wm.ref[:nw][n][:links][a]["length"]
    r = maximum(wm.ref[:nw][n][:resistance][a])

    Δh⁻ = wm.var[:nw][n][:Δh⁻][a]
    q̅ⁿᵉ⁻ = JuMP.upper_bound(wm.var[:nw][n][:q⁻][a])
    rhs⁻ = r * q̅ⁿᵉ⁻^(alpha - 1.0) * wm.var[:nw][n][:q⁻][a]
    con⁻ = JuMP.@constraint(wm.model, inv(L) * Δh⁻ - rhs⁻ <= 0.0)
    wm.con[:nw][n][:directed_potential_loss_ub⁻] = con⁻

    Δh⁺ = wm.var[:nw][n][:Δh⁺][a]
    q̅ⁿᵉ⁺ = JuMP.upper_bound(wm.var[:nw][n][:q⁺][a])
    rhs⁺ = r * q̅ⁿᵉ⁺^(alpha - 1.0) * wm.var[:nw][n][:q⁺][a]
    con⁺ = JuMP.@constraint(wm.model, inv(L) * Δh⁺ - rhs⁺ <= 0.0)
    wm.con[:nw][n][:directed_potential_loss_ub⁺] = con⁺
end

function constraint_link_undirected_flow_ne(wm::GenericWaterModel, a::Int, n::Int=wm.cnw; alpha::Float64=1.852)
    if !haskey(wm.con[:nw][n], :link_undirected_flowⁿᵉ)
        wm.con[:nw][n][:link_undirected_flowⁿᵉ] = Dict{Int, JuMP.ConstraintRef}()
    end

    qⁿᵉ = wm.var[:nw][n][:qⁿᵉ][a]
    q = wm.var[:nw][n][:q][a]
    con⁻ = JuMP.@constraint(wm.model, sum(qⁿᵉ) - q == 0.0)
    wm.con[:nw][n][:link_undirected_flowⁿᵉ] = con
end

function constraint_link_directed_flow_ne(wm::GenericWaterModel, a::Int, n::Int=wm.cnw; alpha::Float64=1.852)
    if !haskey(wm.con[:nw][n], :link_directed_flowⁿᵉ⁻)
        wm.con[:nw][n][:link_directed_flowⁿᵉ⁻] = Dict{Int, JuMP.ConstraintRef}()
        wm.con[:nw][n][:link_directed_flowⁿᵉ⁺] = Dict{Int, JuMP.ConstraintRef}()
    end

    qⁿᵉ⁻ = wm.var[:nw][n][:qⁿᵉ⁻][a]
    q⁻ = wm.var[:nw][n][:q⁻][a]
    con⁻ = JuMP.@constraint(wm.model, sum(qⁿᵉ⁻) - q⁻ == 0.0)
    wm.con[:nw][n][:link_directed_flowⁿᵉ⁻] = con⁻

    qⁿᵉ⁺ = wm.var[:nw][n][:qⁿᵉ⁺][a]
    q⁺ = wm.var[:nw][n][:q⁺][a]
    con⁺ = JuMP.@constraint(wm.model, sum(qⁿᵉ⁺) - q⁺ == 0.0)
    wm.con[:nw][n][:link_directed_flowⁿᵉ⁺] = con⁺
end

function constraint_link_directed_flow(wm::GenericWaterModel, a::Int, n::Int=wm.cnw; alpha::Float64=1.852)
    if !haskey(wm.con[:nw][n], :link_directed_flow)
        wm.con[:nw][n][:link_directed_flow] = Dict{Int, JuMP.ConstraintRef}()
    end

    q = wm.var[:nw][n][:q][a]
    q⁻ = wm.var[:nw][n][:q⁻][a]
    q⁺ = wm.var[:nw][n][:q⁺][a]
    con = JuMP.@constraint(wm.model, (q⁺ - q⁻) - q == 0.0)
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
    xᵈᶦʳ = wm.var[:nw][n][:xᵈᶦʳ]
    out_arcs = filter(a -> i == a.second["node1"], wm.ref[:nw][n][:links])
    out = Array{JuMP.VariableRef}([xᵈᶦʳ[a] for a in keys(out_arcs)])
    in_arcs = filter(a -> i == a.second["node2"], wm.ref[:nw][n][:links])
    in = Array{JuMP.VariableRef}([xᵈᶦʳ[a] for a in keys(in_arcs)])

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
    xᵈᶦʳ = wm.var[:nw][n][:xᵈᶦʳ]
    out_arcs = filter(a -> i == a.second["node1"], wm.ref[:nw][n][:links])
    out = Array{JuMP.VariableRef}([xᵈᶦʳ[a] for a in keys(out_arcs)])
    in_arcs = filter(a -> i == a.second["node2"], wm.ref[:nw][n][:links])
    in = Array{JuMP.VariableRef}([xᵈᶦʳ[a] for a in keys(in_arcs)])

    # Add the sink flow direction constraint.
    con = JuMP.@constraint(wm.model, sum(in) - sum(out) >= 1.0 - length(out))
    wm.con[:nw][n][:directed_sink_flow][i] = con
end
