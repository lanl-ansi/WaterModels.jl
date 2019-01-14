########################################################################
# This file defines commonly-used constraints for water systems models.
########################################################################

function constraint_directed_segmented_flow_conservation(wm::GenericWaterModel, i::Int, n_n::Int)
    # Collect the required variables.
    connections = wm.ref[:nw][n_n][:connection]

    # Initialize the flow sum expression.
    flow_sum = AffExpr(0.0)

    for (a, conn) in filter(a -> i == parse(Int, a.second["node2"]), connections)
        for r in 1:length(wm.ref[:nw][n_n][:resistance][a])
            flow_sum += sum(wm.var[:nw][n_n][:qp][a][:, r])
            flow_sum -= sum(wm.var[:nw][n_n][:qn][a][:, r])
        end
    end

    for (a, conn) in filter(a -> i == parse(Int, a.second["node1"]), connections)
        for r in 1:length(wm.ref[:nw][n_n][:resistance][a])
            flow_sum -= sum(wm.var[:nw][n_n][:qp][a][:, r])
            flow_sum += sum(wm.var[:nw][n_n][:qn][a][:, r])
        end
    end

    if !haskey(wm.con[:nw][n_n], :flow_conservation)
        wm.con[:nw][n_n][:flow_conservation] = Dict{Int, ConstraintRef}()
    end

    # Add the flow conservation constraint.
    demand = wm.ref[:nw][n_n][:junctions][i]["demand"]
    con = @constraint(wm.model, flow_sum == demand)
    wm.con[:nw][n_n][:flow_conservation][i] = con
end

function constraint_directed_flow_conservation(wm::GenericWaterModel, i::Int, n_n::Int)
    # Collect the required variables.
    connections = wm.ref[:nw][n_n][:connection]

    # Initialize the flow sum expression.
    flow_sum = AffExpr(0.0)

    for (a, conn) in filter(a -> i == parse(Int, a.second["node2"]), connections)
        for r in 1:length(wm.ref[:nw][n_n][:resistance][a])
            flow_sum += wm.var[:nw][n_n][:qp][a][r]
            flow_sum -= wm.var[:nw][n_n][:qn][a][r]
        end
    end

    for (a, conn) in filter(a -> i == parse(Int, a.second["node1"]), connections)
        for r in 1:length(wm.ref[:nw][n_n][:resistance][a])
            flow_sum -= wm.var[:nw][n_n][:qp][a][r]
            flow_sum += wm.var[:nw][n_n][:qn][a][r]
        end
    end

    if !haskey(wm.con[:nw][n_n], :flow_conservation)
        wm.con[:nw][n_n][:flow_conservation] = Dict{Int, ConstraintRef}()
    end

    # Add the flow conservation constraint.
    demand = wm.ref[:nw][n_n][:junctions][i]["demand"]
    con = @constraint(wm.model, flow_sum == demand)
    wm.con[:nw][n_n][:flow_conservation][i] = con
end

function constraint_flow_conservation(wm::GenericWaterModel, i::Int, n::Int = wm.cnw)
    # Collect the required variables.
    connections = wm.ref[:nw][n][:connection]

    out_arcs = collect(keys(filter(is_out_node_function(i), connections)))
    out_vars = Array{JuMP.Variable}([wm.var[:nw][n][:q][a] for a in out_arcs])
    in_arcs = collect(keys(filter(is_in_node_function(i), connections)))
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
