######################################################################
# This file defines commonly-used objectives for water systems models.
######################################################################


"""
    objective_wf(wm::AbstractWaterModel)::Nothing

Sets the objective for [Water Flow (WF)](@ref) problem specifications. By
default, only feasibility must be satisfied, i.e., there is no objective.
"""
function objective_wf(wm::AbstractWaterModel)::Nothing
    return JuMP.set_objective_sense(wm.model, JuMP.FEASIBILITY_SENSE)
end


"""
    objective_des(wm::AbstractWaterModel)::JuMP.AffExpr

Sets and returns the objective function for network design (des) problem
specifications. By default, the cost of selecting the discrete pipe resistances
over all design pipes and time indices is minimized. That is, the objective is
```math
\\text{minimize} ~ \\sum_{t \\in \\mathcal{T}}
\\sum_{(i, j) \\in \\mathcal{A}^{\\textrm{des}}_{t}} c_{ijt} z_{ijt},
```
where ``\\mathcal{T}`` is the set of time indices,
``\\mathcal{A}^{\\textrm{des}}_{t}`` is the set of design pipes that are
available at time index ``t``, ``c_{ijt}`` is the cost of installing design
pipe ``(i, j)`` at time index ``t``, and ``z_{ijt}`` is a binary variable
indicating whether (1) or not (0) the design pipe is selected for construction.
"""
function objective_des(wm::AbstractWaterModel)::JuMP.AffExpr
    objective = JuMP.AffExpr(0.0)

    for n in nw_ids(wm)
        for (a, des_pipe) in ref(wm, n, :des_pipe)
            z_des_pipe_term = des_pipe["cost"] * var(wm, n, :z_des_pipe, a)
            JuMP.add_to_expression!(objective, z_des_pipe_term)
        end
    end

    return JuMP.@objective(wm.model, JuMP.MIN_SENSE, objective)
end


"""
    objective_max_demand(wm::AbstractWaterModel)::JuMP.AffExpr
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
    return JuMP.@objective(wm.model, JuMP.MAX_SENSE, objective)
end


"""
    objective_ne(wm::AbstractWaterModel)::JuMP.AffExpr

Sets the objective function for [Network Expansion (ne)](@ref) problem specifications. The
cost of building and operating all network expansion components is minimized.
"""
function objective_ne(wm::AbstractWaterModel)::JuMP.AffExpr
    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(nw_ids(wm)))

    # Find the network IDs over which the objective will be defined.
    if length(network_ids) > 1
        network_ids_status = network_ids[1:end-1]
    else
        network_ids_status = network_ids
    end

    # Initialize the objective expression to zero.
    objective = JuMP.AffExpr(0.0)

    for n in network_ids_status
        # Get the set of network expansion short pipes at time index `n`.
        for (a, ne_short_pipe) in ref(wm, n, :ne_short_pipe)
            # Add the cost of network expansion component `a` at time period `n`.
            term = ne_short_pipe["construction_cost"] * var(wm, n, :z_ne_short_pipe, a)
            JuMP.add_to_expression!(objective, term)
        end

        # Get the set of network expansion pumps at time index `n`.
        for (a, ne_pump) in ref(wm, n, :ne_pump)
            # Add the cost of network expansion component `a` at time period `n`.
            term = ne_pump["construction_cost"] * var(wm, n, :z_ne_pump, a)
            JuMP.add_to_expression!(objective, term)
        end
    end

    # Minimize the total cost of network expansion.
    return JuMP.@objective(wm.model, JuMP.MIN_SENSE, objective)
end


"""
    objective_owf(wm::AbstractWaterModel)::JuMP.AffExpr

Sets the objective for optimal water flow (owf) problem specifications. By
default, minimizes the costs associated with (1) extracting water from each
reservoir at some volumetric flow rate and (2) pumping operations, which
consume power at predefined electricity rates. That is, the objective is
```math
\\text{minimize} ~ \\sum_{t \\in \\mathcal{T}^{\\prime}} \\Delta t
\\left(\\left[\\sum_{(i, j) \\in \\mathcal{P}_{t}} \\pi_{ijt} P_{ijt}\\right] +
\\left[\\sum_{i \\in \\mathcal{R}_{t}} \\mu_{it} q_{it}\\right]\\right),
```
where ``\\mathcal{T}^{\\prime} := \\mathcal{T} \\setminus \\{T\\}`` is the set
of time indices without the last time index, which is reserved only for
computing the final volumes of tanks; ``\\Delta t`` is the time step connecting
time ``t`` with time ``t + 1``; ``\\mathcal{P}_{t}`` is the set of pumps
available for operation at time ``t``; ``\\pi_{ijt}`` is the price of
electricity (cost per unit energy) used for pumping at pump ``(i, j)`` and time
``t``; ``P_{ijt}`` is the power consumption of pump ``(i, j)`` at time ``t``;
``\\mathcal{R}_{t}`` is the set of reservoirs available at time ``t``;
``\\mu_{it}`` is the cost of extracting (and treating) one volumetric unit of
water from reservoir ``i`` at time ``t``; and ``q_{it}`` is the volumetric flow
rate of water extracted from reservoir ``i`` at time ``t``.
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

        for a in ids(wm, n, :ne_pump)
            # Add expansion pump energy costs to the objective.
            JuMP.add_to_expression!(objective, var(wm, n, :c_ne_pump, a))
        end
    end

    # Minimize the cost of network operation.
    return JuMP.@objective(wm.model, JuMP.MIN_SENSE, objective)
end
