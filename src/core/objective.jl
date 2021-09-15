######################################################################
# This file defines commonly-used objectives for water systems models.
######################################################################


"""
    objective_wf(wm::AbstractWaterModel)::Nothing

Sets the objective function for [Water Flow (WF)](@ref) and [Multinetwork Water Flow (MN
WF)](@ref) problem specifications. By default, only feasibility must be satisfied.
"""
function objective_wf(wm::AbstractWaterModel)::Nothing
    return JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end


"""
    objective_des(wm::AbstractWaterModel)::JuMP.AffExpr

Sets the objective function for [Optimal Design (DES)](@ref) problem specifications. The
cost of selecting the discrete pipe resistances over all design pipes is minimized.
"""
function objective_des(wm::AbstractWaterModel)::JuMP.AffExpr
    objective = JuMP.AffExpr(0.0)

    for n in nw_ids(wm)
        for (a, des_pipe) in ref(wm, n, :des_pipe)
            z_des_pipe_term = des_pipe["cost"] * var(wm, n, :z_des_pipe, a)
            JuMP.add_to_expression!(objective, z_des_pipe_term)
        end
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end


"""
    objective_max_demand(wm::AbstractWaterModel)::JuMP.AffExpr

Sets the objective function for [Maximal Demand Delivery (MDD)](@ref) problem specifications.
"""
function objective_max_demand(wm::AbstractWaterModel)::JuMP.AffExpr
    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(nw_ids(wm)))

    # Find the network IDs over which the objective will be defined.
    if length(network_ids) > 1
        network_ids_flow = network_ids[1:end-1]
    else
        network_ids_flow = network_ids
    end

    # Initialize the objective expression to zero.
    objective = JuMP.AffExpr(0.0)

    for n in network_ids_flow
        for (i, demand) in ref(wm, n, :dispatchable_demand)
            # Add the volume delivered at demand `i` and time period `n`.
            coeff = get(demand, "priority", 1.0) * ref(wm, n, :time_step)
            JuMP.add_to_expression!(objective, coeff * var(wm, n, :q_demand, i))
        end
    end

    # Maximize the total amount of water volume delivered.
    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end


"""
    objective_owf(wm::AbstractWaterModel)::JuMP.AffExpr

Sets the objective function for [Optimal Water Flow (OWF)](@ref) and [Multinetwork Optimal
Water Flow (MN OWF)](@ref) problem specifications.
"""
function objective_owf(wm::AbstractWaterModel)::JuMP.AffExpr
    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(nw_ids(wm)))

    # Find the network IDs to define the objective over.
    if length(network_ids) > 1
        network_ids_flow = network_ids[1:end-1]
    else
        network_ids_flow = network_ids
    end

    # Initialize the objective expression to zero.
    objective = JuMP.AffExpr(0.0)

    for n in network_ids_flow
        for (i, reservoir) in ref(wm, n, :reservoir)
            # Add reservoir extraction and treatment costs to the objective.
            @assert haskey(reservoir, "flow_cost") # Ensure a flow cost exists.
            coeff = reservoir["flow_cost"] * ref(wm, n, :time_step)
            JuMP.add_to_expression!(objective, coeff * var(wm, n, :q_reservoir, i))
        end
    end

    for n in network_ids_flow
        for a in ids(wm, n, :pump)
            # Add pump energy costs to the objective.
            JuMP.add_to_expression!(objective, var(wm, n, :c_pump, a))
        end
    end

    # Minimize the cost of network operation.
    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end