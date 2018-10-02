########################################################################
# This file defines commonly-used constraints for water systems models.
########################################################################

function constraint_flow_conservation{T}(wm::GenericWaterModel{T}, i, n::Int = wm.cnw)
    # Add the flow conservation constraints for junction nodes.
    out_arcs = collect(keys(filter((id, pipe) -> pipe["node1"] == i, wm.ref[:nw][n][:pipes])))
    out_vars = Array{JuMP.Variable}([wm.var[:nw][n][:q][a] for a in out_arcs])
    in_arcs = collect(keys(filter((id, pipe) -> pipe["node2"] == i, wm.ref[:nw][n][:pipes])))
    in_vars = Array{JuMP.Variable}([wm.var[:nw][n][:q][a] for a in in_arcs])

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

function constraint_no_good{T}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    yp_solution = getvalue(wm.var[:nw][n][:yp])
    one_vars = Array{JuMP.Variable}([wm.var[:nw][n][:yp][idx[1]] for idx in keys(yp_solution) if yp_solution[idx[1]] > 1.0e-4])
    zero_vars = Array{JuMP.Variable}([wm.var[:nw][n][:yp][idx[1]] for idx in keys(yp_solution) if yp_solution[idx[1]] <= 1.0e-4])
    @constraint(wm.model, sum(zero_vars) - sum(one_vars) >= 1 - length(one_vars))
end
