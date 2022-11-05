######################################################################
# This file defines commonly-used objectives for water systems models.
######################################################################

"""
    objective_wf(wm::AbstractWaterModel)

Sets the objective for [Water Flow (WF)](@ref) problem specifications. By
default, only feasibility must be satisfied, i.e., there is no objective.
"""
function objective_wf(wm::AbstractWaterModel)
    JuMP.set_objective_sense(wm.model, JuMP.FEASIBILITY_SENSE)
end


"""
    objective_des(wm::AbstractWaterModel)

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
function objective_des(wm::AbstractWaterModel)
    objective = JuMP.AffExpr(0.0)

    for (n, nw_ref) in nws(wm)
        for (a, des_pipe) in ref(wm, n, :des_pipe)
            z_des_pipe_term = des_pipe["cost"] * var(wm, n, :z_des_pipe, a)
            JuMP.add_to_expression!(objective, z_des_pipe_term)
        end
    end

    return JuMP.@objective(wm.model, JuMP.MIN_SENSE, objective)
end


"""
    objective_owf(wm::AbstractWaterModel)

Sets the objective for optimal water flow (owf) problem specifications. By
default, minimizes the costs associated with (1) extracting water from each
reservoir at some volumetric flow rate and (2) pumping operations, which
consume power at predefined electricity rates. That is, the objective is
"""
function objective_owf(wm::AbstractWaterModel)
    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(nw_ids(wm)))

    # Initialize the objective expression to zero.
    objective = JuMP.AffExpr(0.0)

    for n in network_ids[1:end-1]
        for (i, reservoir) in ref(wm, n, :reservoir)
            # Add reservoir flow extraction and treatment costs to the objective.
            @assert haskey(reservoir, "flow_cost") # Ensure a flow cost exists.
            coeff = reservoir["flow_cost"] * ref(wm, n, :time_step)
            JuMP.add_to_expression!(objective, coeff * var(wm, n, :q_reservoir, i))
        end
    end

    for n in network_ids[1:end-1]
        for (a, pump) in ref(wm, n, :pump)
            # Add pump energy costs to the objective.
            JuMP.add_to_expression!(objective, var(wm, n, :c_pump, a))
        end
    end

    # Normalize the objective to be more numerically well-behaved.
    pos_coeff = filter(x -> x > 0.0, collect(objective.terms.vals))
    minimum_scalar = length(pos_coeff) > 0 ? minimum(pos_coeff) : 1.0
    objective_scaled = (1.0 / minimum_scalar) * objective

    # Minimize the (numerically scaled) cost required to operate pumps.
    return JuMP.@objective(wm.model, JuMP.MIN_SENSE, objective)
end