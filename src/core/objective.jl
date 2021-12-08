######################################################################
# This file defines commonly-used objectives for water systems models.
######################################################################


"""
    objective_wf(wm::AbstractWaterModel)::Nothing
Sets the objective function for [Water Flow (WF)](@ref) and [Multinetwork Water Flow (MN
WF)](@ref) problem specifications. By default, only feasibility must be satisfied.
"""
function objective_wf(wm::AbstractWaterModel)::Nothing
    #return JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
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
        # Get the set of dispatchable demands at time index `n`.
        dispatchable_demands = ref(wm, n, :dispatchable_demand)

        for (i, demand) in filter(x -> x.second["flow_min"] >= 0.0, dispatchable_demands)
            # Add the volume delivered at demand `i` and time period `n`.
            coeff = get(demand, "priority", 1.0) * ref(wm, n, :time_step)
            JuMP.add_to_expression!(objective, coeff * var(wm, n, :q_demand, i))
        end
    end

    # Maximize the total amount of water volume delivered.
    return JuMP.@objective(wm.model, _MOI.MAX_SENSE, objective)
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
    #for n in 1:6
        for (i, reservoir) in ref(wm, n, :reservoir)
            # Add reservoir extraction and treatment costs to the objective.
            @assert haskey(reservoir, "flow_cost") # Ensure a flow cost exists.
            coeff = reservoir["flow_cost"] * ref(wm, n, :time_step)
            JuMP.add_to_expression!(objective, coeff * var(wm, n, :q_reservoir, i))
        end
    end

    for n in network_ids_flow
    #for n in 1:6
        for a in ids(wm, n, :pump)
            # Add pump energy costs to the objective.
            JuMP.add_to_expression!(objective, var(wm, n, :c_pump, a))
        end
    end
    @show objective 
    # Minimize the cost of network operation.
    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end


"""
    objective_owf(wm::AbstractWaterModel)
Sets the objective function for optimal water flow (owf) problem specifications.
"""
# function objective_owf(wm::AbstractWaterModel)
#     objective = JuMP.AffExpr(0.0)

#     #for (n, nw_ref) in nws(wm)
#     for n in 1:4
#         for (i, reservoir) in ref(wm, n, :reservoir)
#             # Add reservoir flow extraction and treatment costs to the objective.
#             @assert haskey(reservoir, "flow_cost") # Ensure a flow cost exists.
#             coeff = reservoir["flow_cost"] * ref(wm, n, :time_step)
#             JuMP.add_to_expression!(objective, coeff * var(wm, n, :q_reservoir, i))
#         end
#     end

#     #for (n, nw_ref) in nws(wm)
#     for n in 1:4
#         for (a, pump) in ref(wm, n, :pump)
#             # Add pump energy costs to the objective.
#             @assert haskey(pump, "energy_price") # Ensure a price exists.
#             coeff = ref(wm, n, :time_step) * pump["energy_price"]
#             JuMP.add_to_expression!(objective, coeff * var(wm, n, :P_pump, a))
#         end
#     end

#     # Normalize the objective to be more numerically well-behaved.
#     # pos_coeff = filter(x -> x > 0.0, collect(objective.terms.vals))
#     # minimum_scalar = length(pos_coeff) > 0 ? minimum(pos_coeff) : 1.0
#     # objective_scaled = (1.0 / minimum_scalar) * objective

#     # Minimize the (numerically scaled) cost required to operate pumps.
#     #return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective_scaled)
#     @show objective
#     return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
# end


# function objective_owf_decomp(wm::AbstractWaterModel)
#     objective = JuMP.AffExpr(0.0)
#     network_ids = sort(collect(nw_ids(wm)))

#     for n in network_ids[2:end]
#         for (i, reservoir) in ref(wm, n, :reservoir)
#             # Add reservoir flow extraction and treatment costs to the objective.
#             @assert haskey(reservoir, "flow_cost") # Ensure a flow cost exists.
#             coeff = reservoir["flow_cost"] * ref(wm, n, :time_step)
#             JuMP.add_to_expression!(objective, coeff * var(wm, n, :q_reservoir, i))
#         end
#     end

#     for n in network_ids[2:end]
#         for (a, pump) in ref(wm, n, :pump)
#             # Add pump energy costs to the objective.
#             @assert haskey(pump, "energy_price") # Ensure a price exists.
#             coeff = ref(wm, n, :time_step) * pump["energy_price"]
#             JuMP.add_to_expression!(objective, coeff * var(wm, n, :P_pump, a))
#         end
#     end

#     # Normalize the objective to be more numerically well-behaved.
#     pos_coeff = filter(x -> x > 0.0, collect(objective.terms.vals))
#     minimum_scalar = length(pos_coeff) > 0 ? minimum(pos_coeff) : 1.0
#     objective_scaled = (1.0 / minimum_scalar) * objective

#     # Minimize the (numerically scaled) cost required to operate pumps.
#     return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective_scaled)
# end

function objective_owf_decomp_v2(wm::AbstractWaterModel)
    objective = JuMP.AffExpr(0.0)
    network_ids = sort(collect(nw_ids(wm)))

    for n in network_ids[1:end - 1]
        for (i, reservoir) in ref(wm, n, :reservoir)
            # Add reservoir flow extraction and treatment costs to the objective.
            @assert haskey(reservoir, "flow_cost") # Ensure a flow cost exists.
            coeff = reservoir["flow_cost"] * ref(wm, n, :time_step)
            JuMP.add_to_expression!(objective, coeff * var(wm, n, :q_reservoir, i))
        end
    end

    for n in network_ids[1:end - 1]
        for (a, pump) in ref(wm, n, :pump)
            # Add pump energy costs to the objective.
            @assert haskey(pump, "energy_price") # Ensure a price exists.
            coeff = ref(wm, n, :time_step) * pump["energy_price"]
            JuMP.add_to_expression!(objective, coeff * var(wm, n, :P_pump, a))
        end
    end

    # Normalize the objective to be more numerically well-behaved.
    # pos_coeff = filter(x -> x > 0.0, collect(objective.terms.vals))
    # minimum_scalar = length(pos_coeff) > 0 ? minimum(pos_coeff) : 1.0
    # objective_scaled = (1.0 / minimum_scalar) * objective

    # Minimize the (numerically scaled) cost required to operate pumps.
    #return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective_scaled)
    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end


function objective_owf_decomp_v3(wm::AbstractWaterModel)
    objective = JuMP.AffExpr(0.0)
    network_ids = sort(collect(nw_ids(wm)))

    for n in network_ids[1:end - 1]
        for (i, reservoir) in ref(wm, n, :reservoir)
            # Add reservoir flow extraction and treatment costs to the objective.
            @assert haskey(reservoir, "flow_cost") # Ensure a flow cost exists.
            coeff = reservoir["flow_cost"] * ref(wm, n, :time_step)
            JuMP.add_to_expression!(objective, coeff * var(wm, n, :q_reservoir, i))
        end
    end

    for n in network_ids[1:end - 1]
        for (a, pump) in ref(wm, n, :pump)
            # Add pump energy costs to the objective.
            @assert haskey(pump, "energy_price") # Ensure a price exists.
            coeff = ref(wm, n, :time_step) * pump["energy_price"]
            JuMP.add_to_expression!(objective, coeff * var(wm, n, :P_pump, a))
        end
    end

    # Add dual term to the initial tank position
    for (t, tank) in ref(wm, network_ids[1], :tank)
        @assert haskey(tank, "initial_tank_h")
        @assert haskey(tank, "initial_tank_dual")

        lambda = tank["initial_tank_dual"]
        h0     = tank["initial_tank_h"]
        JuMP.add_to_expression!(objective, lambda * (var(wm, network_ids[1], :h, t) - h0) )
    end

    #Add dual term to the final tank position
    for (t, tank) in ref(wm, netowrk_ids[end], :tank)
        @assert haskey(tank, "final_tank_h")
        @assert haskey(tank, "final_tank_dual")

        lambda = tank["final_tank_dual"]
        he     = tank["final_tank_h"]
        JuMP.add_to_expression!(objective, lambda * (var(wm, netowrk_ids[end], :h, t) - he) )
    end

    # # Normalize the objective to be more numerically well-behaved.
    # pos_coeff = filter(x -> x > 0.0, collect(objective.terms.vals))
    # minimum_scalar = length(pos_coeff) > 0 ? minimum(pos_coeff) : 1.0
    # objective_scaled = (1.0 / minimum_scalar) * objective

    # Minimize the (numerically scaled) cost required to operate pumps.
    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end



function objective_geo_owf(wm::AbstractWaterModel)
    objective = JuMP.AffExpr(0.0)
    network_ids = sort(collect(nw_ids(wm)))

    for n in network_ids[1:end - 1]
        for (i, reservoir) in ref(wm, n, :reservoir)
            # Add reservoir flow extraction and treatment costs to the objective.
            @assert haskey(reservoir, "flow_cost") # Ensure a flow cost exists.
            coeff = reservoir["flow_cost"] * ref(wm, n, :time_step)
            JuMP.add_to_expression!(objective, coeff * var(wm, n, :q_reservoir, i))
        end
    end

    for n in network_ids[1:end - 1]
        for (a, pump) in ref(wm, n, :pump)
            # Add pump energy costs to the objective.
            @assert haskey(pump, "energy_price") # Ensure a price exists.
            coeff = ref(wm, n, :time_step) * pump["energy_price"]
            JuMP.add_to_expression!(objective, coeff * var(wm, n, :P_pump, a))
        end
    end

    # for n in network_ids[1:end-1]
    #     cmn_comp = values(wm.data["nw"][string(n)]["common"])
    #
    #     for comp_idx in cmn_comp
    #         comp = ref(wm, n, :node, parse(Int, comp_idx))
    #         @assert haskey(comp, "head_pressure_val")
    #         @assert haskey(comp, "head_pressure_dual")
    #
    #         h_     = comp["head_pressure_val"]
    #         lambda = comp["head_pressure_dual"]
    # 
    #         JuMP.add_to_expression!(objective, lambda* (var(wm, n, :h, parse(Int, comp_idx)) - h_))
    #     end
    # end


    # Normalize the objective to be more numerically well-behaved.
    # pos_coeff = filter(x -> x > 0.0, collect(objective.terms.vals))
    # minimum_scalar = length(pos_coeff) > 0 ? minimum(pos_coeff) : 1.0
    # objective_scaled = (1.0 / minimum_scalar) * objective

    # Minimize the (numerically scaled) cost required to operate pumps.
    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end