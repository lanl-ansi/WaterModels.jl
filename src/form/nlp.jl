# Define NLP implementations of water distribution models.

export NLPWaterModel, StandardNLPForm

""
@compat abstract type AbstractNLPForm <: AbstractWaterFormulation end

""
@compat abstract type StandardNLPForm <: AbstractNLPForm end

"Default (nonconvex) NLP model."
const NLPWaterModel = GenericWaterModel{StandardNLPForm}

"Default NLP constructor."
NLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardNLPForm; kwargs...)

"Create variables associated with the head for the nonconvex NLP problem."
function variable_head{T <: StandardNLPForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    variable_head_common(wm, n)
end

"Non-convex Darcy-Weisbach constraint with unknown direction."
function constraint_dw_unknown_direction{T <: StandardNLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Add a nonlinear constraint for the head loss.
    @NLconstraint(wm.model, h_i - h_j == lambda * q * abs(q))
end

"Non-convex Darcy-Weisbach constraint with known direction."
function constraint_dw_known_direction{T <: StandardNLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a non-convex quadratic constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) == lambda * q^2)
end

"Non-convex Hazen-Williams constraint for flow with unknown direction."
function constraint_hw_unknown_direction{T <: StandardNLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Add a non-convex constraint for the head loss.
    @NLconstraint(wm.model, h_i - h_j == lambda * q * (q^2)^0.426)
end

"Non-convex Hazen-Williams constraint for flow with known direction."
function constraint_hw_known_direction{T <: StandardNLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a non-convex constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) == lambda * (q^2)^0.926)
end
