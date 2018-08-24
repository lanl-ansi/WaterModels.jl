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

    if haskey(wm.ref[:nw][n][:ne_pipe], a)
        # Get variables associated with the discrete diameter choices.
        diameter_vars = wm.var[:nw][n][:diameter][a]
        diameters = [key[1] for key in keys(diameter_vars)]

        # Get the corresponding diameter measurements (meters).
        pipe = wm.ref[:nw][n][:ne_pipe][a]
        lambdas = [calc_friction_factor_hw_ne(pipe, d) for d in diameters]
        aff = AffExpr(diameter_vars[:], lambdas, 0.0)

        aux = @variable(wm.model)
        @constraint(wm.model, aux == aff)
        @NLconstraint(wm.model, h_i - h_j == aux * q * (q^2)^0.426)

        # Add a constraint that says at most one diameter may be selected.
        @constraint(wm.model, sum(diameter_vars) == 1)
    else
        # Add a non-convex constraint for the head loss.
        lambda = calc_friction_factor_hw(wm.ref[:nw][n][:pipes][a])
        @NLconstraint(wm.model, h_i - h_j == lambda * q * (q^2)^0.426)
    end
end
