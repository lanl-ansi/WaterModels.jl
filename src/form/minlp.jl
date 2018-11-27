# Defines MINLP implementations of water distribution models.

export MINLPWaterModel, StandardMINLPForm

abstract type AbstractMINLPForm <: AbstractWaterFormulation end
abstract type StandardMINLPForm <: AbstractMINLPForm end

"Default (nonconvex) MINLP model."
const MINLPWaterModel = GenericWaterModel{StandardMINLPForm}

"Default MINLP constructor."
MINLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMINLPForm; kwargs...)

function constraint_flow_direction(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: StandardMINLPForm
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)
end

#"Non-convex Darcy-Weisbach constraint with unknown direction."
#function constraint_dw_unknown_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMINLPForm
#    # Collect variables and parameters needed for the constraint.
#    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)
#
#    # Add a nonlinear constraint for the head loss.
#    @NLconstraint(wm.model, h_i - h_j == lambda * q * abs(q))
#end
#
#"Non-convex Hazen-Williams constraint for flow with unknown direction."
#function constraint_hw_unknown_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMINLPForm
#    # Collect variables and parameters needed for the constraint.
#    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)
#
#    # Add a non-convex constraint for the head loss.
#    @NLconstraint(wm.model, h_i - h_j == lambda * q * (q^2)^0.426)
#end
#
