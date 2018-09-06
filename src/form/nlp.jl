# Defines NLP implementations of water distribution models.

export NLPWaterModel, StandardNLPForm

""
@compat abstract type AbstractNLPForm <: AbstractWaterFormulation end

""
@compat abstract type StandardNLPForm <: AbstractNLPForm end

"Default (nonconvex) NLP model."
const NLPWaterModel = GenericWaterModel{StandardNLPForm}

"Default NLP constructor."
NLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardNLPForm; kwargs...)

"Non-convex Darcy-Weisbach constraint with unknown direction."
function constraint_dw_unknown_direction{T <: StandardNLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Add a nonlinear constraint for the head loss.
    @NLconstraint(wm.model, h_i - h_j == lambda * q * abs(q))
end

"Non-convex Hazen-Williams constraint for flow with unknown direction."
function constraint_hw_unknown_direction{T <: StandardNLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Add a non-convex constraint for the head loss.
    @NLconstraint(wm.model, h_i - h_j == lambda * q * (q^2)^0.426)
end

"Non-convex Hazen-Williams constraint for flow with unknown direction."
function constraint_hw_unknown_direction_ne{T <: StandardNLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

    # Add constraints required to define gamma.
    constraint_define_gamma_hw_ne(wm, a, n)

    # Define an auxiliary variable for the sum of the gamma variables.
    gamma_sum = @variable(wm.model, basename = "gamma_sum_$(n)_$(a)", start = 0.0)
    @constraint(wm.model, gamma_sum == sum(wm.var[:nw][n][:gamma][a]))

    # Add a non-convex constraint for the head loss.
    @NLconstraint(wm.model, gamma_sum == q * (q^2)^0.426)
end
