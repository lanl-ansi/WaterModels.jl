# Define MICP (mixed-integer convex program) implementations of water distribution models.

export MICPWaterModel, StandardMICPForm

abstract type AbstractMICPForm <: AbstractWaterFormulation end
abstract type StandardMICPForm <: AbstractMICPForm end

"The default MICP (mixed-integer convex program) model is a relaxation of the non-convex MINLP model."
const MICPWaterModel = GenericWaterModel{StandardMICPForm}

"Default MICP constructor."
MICPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMICPForm; kwargs...)

"Convex (relaxed) Darcy-Weisbach constraint for flow with unknown direction."
function constraint_dw_unknown_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: AbstractMICPForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Add constraints to restrict the direction.
    constraint_restrict_direction(wm, a, n)

    # Add constraints required to define gamma.
    constraint_define_gamma(wm, a, n)

    # Add a convex quadratic constraint for the head loss.
    gamma = wm.var[:nw][n][:gamma][a]
    @NLconstraint(wm.model, gamma >= lambda * q^2)
end

"Convex (relaxed) Hazen-Williams constraint for flow with unknown direction."
function constraint_hw_unknown_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: AbstractMICPForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Add constraints to restrict the direction.
    constraint_restrict_direction(wm, a, n)

    # Add constraints required to define gamma.
    constraint_define_gamma(wm, a, n)

    # Add a nonlinear constraint for the head loss.
    gamma = wm.var[:nw][n][:gamma][a]
    @NLconstraint(wm.model, gamma >= lambda * (q^2)^0.926)
end

"Convex (relaxed) Hazen-Williams constraint for flow with unknown direction."
function constraint_hw_unknown_direction_ne(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: AbstractMICPForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

    # Add constraints to restrict the direction.
    constraint_restrict_direction(wm, a, n)

    # Add constraints required to define gamma.
    constraint_define_gamma_hw_ne(wm, a, n)

    # Define an auxiliary variable for the sum of the gamma variables.
    gamma_sum = wm.var[:nw][n][:gamma_sum][a]
    @constraint(wm.model, gamma_sum == sum(wm.var[:nw][n][:gamma][a]))

    # Add a nonlinear constraint for the head loss.
    @NLconstraint(wm.model, gamma_sum >= (q^2)^0.926)
end

"Convex (relaxed) Darcy-Weisbach constraint for flow with unknown direction."
function constraint_dw_unknown_direction_ne(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: AbstractMICPForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

    # Add constraints to restrict the direction.
    constraint_restrict_direction(wm, a, n)

    # Add constraints required to define gamma.
    constraint_define_gamma_dw_ne(wm, a, n)

    # Define an auxiliary variable for the sum of the gamma variables.
    gamma_sum = wm.var[:nw][n][:gamma_sum][a]
    @constraint(wm.model, gamma_sum == sum(wm.var[:nw][n][:gamma][a]))

    # Add a nonlinear constraint for the head loss.
    @NLconstraint(wm.model, gamma_sum >= q^2)
end


"Convex (relaxed) Darcy-Weisbach constraint for flow with known direction."
function constraint_dw_known_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: AbstractMICPForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a nonlinear constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) >= lambda * q^2)
end

"Convex (relaxed) Hazen-Williams constraint for flow with known direction."
function constraint_hw_known_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: AbstractMICPForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a nonlinear constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) >= lambda * (q^2)^0.926)
end
