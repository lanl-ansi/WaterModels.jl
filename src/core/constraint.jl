########################################################################
# This file defines commonly-used constraints for water systems models.
########################################################################

function constraint_flow_conservation{T}(wm::GenericWaterModel{T}, i, n::Int = wm.cnw)
    # Collect the required variables.
    connections = wm.ref[:nw][n][:connection]
    out_arcs = collect(keys(filter((id, c) -> c["node1"] == i, connections)))
    out_vars = Array{JuMP.Variable}([wm.var[:nw][n][:q][a] for a in out_arcs])
    in_arcs = collect(keys(filter((id, c) -> c["node2"] == i, connections)))
    in_vars = Array{JuMP.Variable}([wm.var[:nw][n][:q][a] for a in in_arcs])

    # Add the flow conservation constraints for junction nodes.
    if haskey(wm.ref[:nw][n][:junctions], i)
        demand = wm.ref[:nw][n][:junctions][i]["demand"]
        @constraint(wm.model, sum(in_vars) - sum(out_vars) == demand)
    elseif haskey(wm.ref[:nw][n][:reservoirs], i)
        junctions = values(wm.ref[:nw][n][:junctions])
        sum_demand = sum(junction["demand"] for junction in junctions)
        @constraint(wm.model, sum(out_vars) - sum(in_vars) >= 0.0)
        @constraint(wm.model, sum(out_vars) - sum(in_vars) <= sum_demand)
    end
end
