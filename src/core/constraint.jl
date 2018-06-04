########################################################################
# This file defines commonly-used constraints for water systems models.
########################################################################

function compute_lambda(ref, id)
    # Diameter assumes original units of millimeters.
    diameter = ref[:pipes][id]["diameter"] / 1000.0

    # Roughness assumes original units of millimeters.
    roughness = ref[:pipes][id]["roughness"] / 1000.0

    # Length assumes original units of meters.
    length = ref[:pipes][id]["length"]

    # Use standard gravitational acceleration on earth.
    g = 9.80665

    # Use the Prandtl-Kármán friction factor.
    f_s = 0.25 / log((roughness / diameter) / 3.71)^2

    # Return lambda.
    return (8.0 * length) / (pi^2 * g * diameter^5) * f_s
end

function compute_hazen_williams_lambda(ref, id)
    # Diameter assumes original units of millimeters.
    diameter = ref[:pipes][id]["diameter"] / 1000.0

    # Roughness assumes no units (?).
    roughness = ref[:pipes][id]["roughness"] / 1000.0

    # Length assumes original units of meters.
    length = ref[:pipes][id]["length"]

    # Return lambda.
    return (10.67 * length) / (roughness^1.852 * diameter^4.87)
end

function constraint_flow_conservation{T}(wm::GenericWaterModel{T}, i, n::Int = wm.cnw)
    # Add the flow conservation constraints for junction nodes.
    out_arcs = collect(keys(filter((id, pipe) -> pipe["node1"] == i, wm.ref[:nw][n][:pipes])))
    out_vars = Array{JuMP.Variable}([wm.var[:nw][n][:q][a] for a in out_arcs])
    in_arcs = collect(keys(filter((id, pipe) -> pipe["node2"] == i, wm.ref[:nw][n][:pipes])))
    in_vars = Array{JuMP.Variable}([wm.var[:nw][n][:q][a] for a in in_arcs])

    # Demands assume original units of liters per second.
    if haskey(wm.ref[:nw][n][:junctions], i)
        demand = 0.001 * wm.ref[:nw][n][:junctions][i]["demand"]
        @constraint(wm.model, sum(in_vars) - sum(out_vars) == demand)
    elseif haskey(wm.ref[:nw][n][:reservoirs], i)
        all_demands = [junction["demand"] for junction in values(wm.ref[:nw][n][:junctions])]
        sum_demand = 0.001 * sum(all_demands)
        #@constraint(wm.model, sum(out_vars) - sum(in_vars) >= 0.0)
        @constraint(wm.model, sum(out_vars) - sum(in_vars) >= sum_demand)
    end
end

function constraint_potential_flow_coupling{T}(wm::GenericWaterModel{T}, i, n::Int = wm.cnw)
    # Get the outgoing arcs of the node.
    arcs_from = collect(keys(filter((id, pipe) -> pipe["node1"] == i, wm.ref[:nw][n][:pipes])))

    for a in arcs_from
        # Collect variables needed for the constraint.
        gamma = wm.var[:nw][n][:gamma][a]
        q = wm.var[:nw][n][:q][a]

        #lambda = compute_lambda(wm.ref[:nw][n], a)
        lambda = compute_hazen_williams_lambda(wm.ref[:nw][n], a)
        @constraint(wm.model, gamma >= lambda * q * q)
    end
end

function constraint_define_gamma{T}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Get source and target nodes corresponding to the arc.
    i = wm.ref[:nw][n][:pipes][a]["node1"]
    j = wm.ref[:nw][n][:pipes][a]["node2"]

    # Collect variables needed for the constraint.
    gamma = wm.var[:nw][n][:gamma][a]
    y_p = wm.var[:nw][n][:yp][a]
    y_n = wm.var[:nw][n][:yn][a]
    h_i = wm.var[:nw][n][:h][i]
    h_j = wm.var[:nw][n][:h][j]

    # Get additional data related to the variables.
    h_i_lb = getlowerbound(h_i)
    h_i_ub = getupperbound(h_i)
    h_j_lb = getlowerbound(h_j)
    h_j_ub = getupperbound(h_j)

    # Add the required constraints.
    @constraint(wm.model, h_j - h_i + (h_i_lb - h_j_ub) <= gamma)
    @constraint(wm.model, h_i - h_j + (h_i_ub - h_j_lb) <= gamma)
    @constraint(wm.model, h_j - h_i + (h_i_ub - h_j_lb) >= gamma)
    @constraint(wm.model, h_i - h_j + (h_i_lb - h_j_ub) >= gamma)
end

function constraint_bidirectional_flow{T}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Get source and target nodes corresponding to the arc.
    i = wm.ref[:nw][n][:pipes][a]["node1"]
    j = wm.ref[:nw][n][:pipes][a]["node2"]

    # Collect variables needed for the constraint.
    q = wm.var[:nw][n][:q][a]
    y_p = wm.var[:nw][n][:yp][a]
    y_n = wm.var[:nw][n][:yn][a]
    h_i = wm.var[:nw][n][:h][i]
    h_j = wm.var[:nw][n][:h][j]

    # Get additional data related to the variables.
    h_i_lb = getlowerbound(h_i)
    h_i_ub = getupperbound(h_i)
    h_j_lb = getlowerbound(h_j)
    h_j_ub = getupperbound(h_j)

    # Add the first set of constraints.
    all_demands = [junction["demand"] for junction in values(wm.ref[:nw][n][:junctions])]
    sum_demand = 0.001 * sum(all_demands)
    @constraint(wm.model, (y_p - 1) * sum_demand <= q)
    @constraint(wm.model, (1 - y_n) * sum_demand >= q)

    # Add the second set of constraints.
    @constraint(wm.model, (1 - y_p) * (h_i_lb - h_j_ub) <= h_i - h_j)
    @constraint(wm.model, (1 - y_n) * (h_i_ub - h_j_lb) >= h_i - h_j)

    # Add the third set of constraints.
    @constraint(wm.model, y_p + y_n == 1)
end

# Flow conservation constraint at a junction.
#"`p[f_idx]^2 + q[f_idx]^2 <= rate_a^2`"
#function constraint_flow_conservation(wm::GenericWaterModel, id::String, f_idx, rate_a)
#    p_fr = pm.var[:nw][n][:p][f_idx]
#    q_fr = pm.var[:nw][n][:q][f_idx]
#    @constraint(pm.model, p_fr^2 + q_fr^2 <= rate_a^2)
#end


#" Constraint that states a flow direction must be chosen "
#function constraint_flow_direction_choice{T}(wm::GenericWaterModel{T}, i)
#    yp = wm.var[:yp][i]
#    yn = wm.var[:yn][i]
#
#    c = @constraint(wm.model, yp + yn == 1)
#
#    if !haskey(wm.constraint, :flow_direction_choice)
#        wm.constraint[:flow_direction_choice] = Dict{Int,ConstraintRef}()
#    end
#    wm.constraint[:flow_direction_choice][i] = c
#end
#
#
#
#" constraints on head drop across pipes "
#function constraint_on_off_head_drop{T}(wm::GenericWaterModel{T}, pipe_idx)
#
#    pipe = wm.ref[:connection][pipe_idx]
#    i_junction_idx = pipe["f_junction"]
#    j_junction_idx = pipe["t_junction"]
#
#    yp = wm.var[:yp][pipe_idx]
#    yn = wm.var[:yn][pipe_idx]
#
#    hi = wm.var[:h][i_junction_idx]
#    hj = wm.var[:h][j_junction_idx]
#
#    hd_min = pipe["hd_min"]
#    hd_max = pipe["hd_max"]
#
#    c1 = @constraint(wm.model, (1-yp) * hd_min <= hi - hj)
#    c2 = @constraint(wm.model, hi - hj <= (1-yn)* hd_max)
#
#    if !haskey(wm.constraint, :on_off_head_drop1)
#        wm.constraint[:on_off_head_drop1] = Dict{Int,ConstraintRef}()
#        wm.constraint[:on_off_head_drop2] = Dict{Int,ConstraintRef}()
#    end
#    wm.constraint[:on_off_head_drop1][pipe_idx] = c1
#    wm.constraint[:on_off_head_drop2][pipe_idx] = c2
#end
#
#" constraints on head due to elevation "
#function constraint_elevation_bound_regular{T}(wm::GenericWaterModel{T}, i)
#    junction = wm.ref[:junction][i]
#
#    h = wm.var[:h][i]
#
#
#    c = @constraint(wm.model, h >= junction["head"] )
#
#    if !haskey(wm.constraint, :elevation_bound_reg)
#        wm.constraint[:elevation_bound_reg] = Dict{Int,ConstraintRef}()
#    end
#    wm.constraint[:elevation_bound_reg][i] = c
#end
#
#" constraints on head due to elevation "
#function constraint_elevation_bound_source{T}(wm::GenericWaterModel{T}, i)
#    junction = wm.ref[:junction][i]
#
#    h = wm.var[:h][i]
#
#    c = @constraint(wm.model, h == junction["head"] )
#
#    if !haskey(wm.constraint, :elevation_bound_src)
#        wm.constraint[:elevation_bound_src] = Dict{Int,ConstraintRef}()
#    end
#    wm.constraint[:elevation_bound_src][i] = c
#end
#
#
#" constraints on flow across pipes "
#function constraint_on_off_pipe_flow_direction{T}(wm::GenericWaterModel{T}, pipe_idx)
#    pipe = wm.ref[:connection][pipe_idx]
#
#    i_junction_idx = pipe["f_junction"]
#    j_junction_idx = pipe["t_junction"]
#
#    yp = wm.var[:yp][pipe_idx]
#    yn = wm.var[:yn][pipe_idx]
#    f  = wm.var[:f][pipe_idx]
#
#    max_flow = wm.ref[:max_flow]
#    hd_max = pipe["hd_max"]
#    hd_min = pipe["hd_min"]
#    w = pipe["resistance"]
#
#    c1 = @constraint(wm.model, -(1-yp)*min(max_flow, sqrt(w*max(hd_max, abs(hd_min)))) <= f)
#    c2 = @constraint(wm.model, f <= (1-yn)*min(max_flow, sqrt(w*max(hd_max, abs(hd_min)))))
#
#    if !haskey(wm.constraint, :on_off_pipe_flow_direction1)
#        wm.constraint[:on_off_pipe_flow_direction1] = Dict{Int,ConstraintRef}()
#        wm.constraint[:on_off_pipe_flow_direction2] = Dict{Int,ConstraintRef}()
#    end
#    wm.constraint[:on_off_pipe_flow_direction1][pipe_idx] = c1
#    wm.constraint[:on_off_pipe_flow_direction2][pipe_idx] = c2
#end
#
#" constraints on flow across pipes due to pipe diameter"
#function constraint_on_off_pipe_flow_direction_diameter{T}(wm::GenericWaterModel{T}, pipe_idx)
#    pipe = wm.ref[:connection][pipe_idx]
#
#    i_junction_idx = pipe["f_junction"]
#    j_junction_idx = pipe["t_junction"]
#    D = pipe["diameter"]
#
#    yp = wm.var[:yp][pipe_idx]
#    yn = wm.var[:yn][pipe_idx]
#    f  = wm.var[:f][pipe_idx]
#
#    max_flow = (pi/4)*(D^2)
#    hd_max = pipe["hd_max"]
#    hd_min = pipe["hd_min"]
#
#
#    c1 = @constraint(wm.model, -(1-yp)*max_flow <= f)
#    c2 = @constraint(wm.model, f <= (1-yn)*max_flow)
#
#    if !haskey(wm.constraint, :on_off_pipe_flow_direction1)
#        wm.constraint[:on_off_pipe_flow_direction_diameter1] = Dict{Int,ConstraintRef}()
#        wm.constraint[:on_off_pipe_flow_direction_diameter2] = Dict{Int,ConstraintRef}()
#    end
#    wm.constraint[:on_off_pipe_flow_direction_diameter1][pipe_idx] = c1
#    wm.constraint[:on_off_pipe_flow_direction_diameter2][pipe_idx] = c2
#end
#
#" standard flow balance equation where demand is fixed "
#function constraint_junction_flow_balance{T}(wm::GenericWaterModel{T}, i)
#    junction = wm.ref[:junction][i]
#    junction_branches = wm.ref[:junction_connections][i]
#
#    f_branches = collect(keys(filter( (a, connection) -> connection["f_junction"] == i, wm.ref[:connection])))
#    t_branches = collect(keys(filter( (a, connection) -> connection["t_junction"] == i, wm.ref[:connection])))
#
#    h = wm.var[:h]
#    f = wm.var[:f]
#
#    c = @constraint(wm.model, junction["demand"] == sum(f[a] for a in t_branches) - sum(f[a] for a in f_branches) )
#
#    if !haskey(wm.constraint, :junction_flow_balance)
#        wm.constraint[:junction_flow_balance] = Dict{Int,ConstraintRef}()
#    end
#    wm.constraint[:junction_flow_balance][i] = c
#end
#
#
#" Make sure there is at least one direction set to take flow away from a junction (typically used on source nodes) "
#function constraint_source_flow{T}(wm::GenericWaterModel{T}, i)
#    f_branches = collect(keys(filter( (a,connection) -> connection["f_junction"] == i, wm.ref[:connection])))
#    t_branches = collect(keys(filter( (a,connection) -> connection["t_junction"] == i, wm.ref[:connection])))
#
#    yp = wm.var[:yp]
#    yn = wm.var[:yn]
#
#    if length(t_branches)!=0 && length(f_branches)!=0
#      c = @constraint(wm.model, sum(yp[a] for a in f_branches) + sum(yn[a] for a in t_branches) >= 1)
#    end
#
#    if length(t_branches)==0 && length(f_branches)!=0
#      c = @constraint(wm.model, sum(yp[a] for a in f_branches)  >= 1)
#    end
#
#    if length(t_branches)!=0 && length(f_branches)==0
#      c = @constraint(wm.model, sum(yn[a] for a in t_branches) >= 1)
#    end
#
#
#    if !haskey(wm.constraint, :source_flow)
#        wm.constraint[:source_flow] = Dict{Int,ConstraintRef}()
#    end
#    wm.constraint[:source_flow][i] = c
#end
#
#
#" Make sure there is at least one direction set to take flow to a junction (typically used on sink nodes) "
#function constraint_sink_flow{T}(wm::GenericWaterModel{T}, i)
#    f_branches = collect(keys(filter( (a,connection) -> connection["f_junction"] == i, wm.ref[:connection])))
#    t_branches = collect(keys(filter( (a,connection) -> connection["t_junction"] == i, wm.ref[:connection])))
#
#    yp = wm.var[:yp]
#    yn = wm.var[:yn]
#
#    c = @constraint(wm.model, sum(yn[a] for a in f_branches) + sum(yp[a] for a in t_branches) >= 1)
#    if !haskey(wm.constraint, :sink_flow)
#        wm.constraint[:sink_flow] = Dict{Int,ConstraintRef}()
#    end
#    wm.constraint[:sink_flow][i] = c
#end
#
#
#" This constraint is intended to ensure that flow is on direction through a node with degree 2 and no production or consumption "
#function constraint_conserve_flow{T}(wm::GenericWaterModel{T}, idx)
#    first = nothing
#    last = nothing
#
#    for i in wm.ref[:junction_connections][idx]
#        connection = wm.ref[:connection][i]
#        if connection["f_junction"] == idx
#            other = connection["t_junction"]
#        else
#            other = connection["f_junction"]
#        end
#
#        if first == nothing
#            first = other
#        elseif first != other
#            if last != nothing && last != other
#                error(string("Error: adding a degree 2 constraint to a node with degree > 2: Junction ", idx))
#            end
#            last = other
#        end
#    end
#
#    yp_first = filter(i -> wm.ref[:connection][i]["f_junction"] == first, wm.ref[:junction_connections][idx])
#    yn_first = filter(i -> wm.ref[:connection][i]["t_junction"] == first, wm.ref[:junction_connections][idx])
#    yp_last  = filter(i -> wm.ref[:connection][i]["t_junction"] == last,  wm.ref[:junction_connections][idx])
#    yn_last  = filter(i -> wm.ref[:connection][i]["f_junction"] == last,  wm.ref[:junction_connections][idx])
#
#    yp = wm.var[:yp]
#    yn = wm.var[:yn]
#
#    c1 = nothing
#    c2 = nothing
#    c3 = nothing
#    c4 = nothing
#    if length(yn_first) > 0 && length(yp_last) > 0
#        for i1 in yn_first
#            for i2 in yp_last
#                c1 = @constraint(wm.model, yn[i1]  == yp[i2])
#                c2 = @constraint(wm.model, yp[i1]  == yn[i2])
#                c3 = @constraint(wm.model, yn[i1] + yn[i2] == 1)
#                c4 = @constraint(wm.model, yp[i1] + yp[i2] == 1)
#            end
#        end
#    end
#
#
#   if length(yn_first) > 0 && length(yn_last) > 0
#        for i1 in yn_first
#            for i2 in yn_last
#                c1 = @constraint(wm.model, yn[i1] == yn[i2])
#                c2 = @constraint(wm.model, yp[i1] == yp[i2])
#                c3 = @constraint(wm.model, yn[i1] + yp[i2] == 1)
#                c4 = @constraint(wm.model, yp[i1] + yn[i2] == 1)
#            end
#        end
#    end
#
#    if length(yp_first) > 0 && length(yp_last) > 0
#        for i1 in yp_first
#            for i2 in yp_last
#                c1 = @constraint(wm.model, yp[i1]  == yp[i2])
#                c2 = @constraint(wm.model, yn[i1]  == yn[i2])
#                c3 = @constraint(wm.model, yp[i1] + yn[i2] == 1)
#                c4 = @constraint(wm.model, yn[i1] + yp[i2] == 1)
#            end
#        end
#    end
#
#    if length(yp_first) > 0 && length(yn_last) > 0
#        for i1 in yp_first
#            for i2 in yn_last
#                c1 = @constraint(wm.model, yp[i1] == yn[i2])
#                c2 = @constraint(wm.model, yn[i1] == yp[i2])
#                c3 = @constraint(wm.model, yp[i1] + yp[i2] == 1)
#                c4 = @constraint(wm.model, yn[i1] + yn[i2] == 1)
#            end
#        end
#    end
#
#    if !haskey(wm.constraint, :conserve_flow1)
#        wm.constraint[:conserve_flow1] = Dict{Int,ConstraintRef}()
#        wm.constraint[:conserve_flow2] = Dict{Int,ConstraintRef}()
#        wm.constraint[:conserve_flow3] = Dict{Int,ConstraintRef}()
#        wm.constraint[:conserve_flow4] = Dict{Int,ConstraintRef}()
#    end
#
#    wm.constraint[:conserve_flow1][idx] = c1
#    wm.constraint[:conserve_flow2][idx] = c2
#    wm.constraint[:conserve_flow3][idx] = c3
#    wm.constraint[:conserve_flow4][idx] = c4
#end
#
#
#" ensures that parallel lines have flow in the same direction "
#function constraint_parallel_flow{T}(wm::GenericWaterModel{T}, idx)
#    connection = wm.ref[:connection][idx]
#    i = min(connection["f_junction"], connection["t_junction"])
#    j = max(connection["f_junction"], connection["t_junction"])
#
#    f_connections = filter(i -> wm.ref[:connection][i]["f_junction"] == connection["f_junction"], wm.ref[:parallel_connections][(i,j)])
#    t_connections = filter(i -> wm.ref[:connection][i]["f_junction"] != connection["f_junction"], wm.ref[:parallel_connections][(i,j)])
#
#    yp = wm.var[:yp]
#    yn = wm.var[:yn]
#
#    if length(wm.ref[:parallel_connections][(i,j)]) <= 1
#        return nothing
#    end
#
#    c = @constraint(wm.model, sum(yp[i] for i in f_connections) + sum(yn[i] for i in t_connections) == yp[idx] * length(wm.ref[:parallel_connections][(i,j)]))
#    if !haskey(wm.constraint, :parallel_flow)
#        wm.constraint[:parallel_flow] = Dict{Int,ConstraintRef}()
#    end
#    wm.constraint[:parallel_flow][idx] = c
#end
