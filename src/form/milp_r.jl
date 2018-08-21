# Define MILP implementations of water distribution models.

using PiecewiseLinearOpt

export MILPWaterModel, RelaxedMILPWaterModel, StandardMILPForm, RelaxedMILPForm

"AbstractMILPForm is derived from AbstractMINLPForm"
@compat abstract type AbstractMILPForm <: AbstractMINLPForm end

"StandardMILPForm is derived from StandardMINLPForm"
@compat abstract type StandardMILPForm <: StandardMINLPForm end

"RelaxedMILP is derived from MICPForm"
@compat abstract type RelaxedMILPForm <: StandardMICPForm end

"The default MILP model is assumed to meet all constraints with equality."
const MILPWaterModel = GenericWaterModel{StandardMILPForm}

"The relaxed MILP model relaxes the head loss constraint to an inequality."
const RelaxedMILPWaterModel = GenericWaterModel{RelaxedMILPForm}

"Default MILP constructor (assumes equality constraints)."
MILPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMILPForm; kwargs...)

"Relaxed MILP constructor."
RelaxedMILPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, RelaxedMILPForm; kwargs...)

"Return values that approximate the Hazen-Williams head loss constraint."
function construct_hw_separators(q::JuMP.Variable, lambda::Float64, n::Int = 3)
    q_points = linspace(getlowerbound(q), getupperbound(q), n)
    f_evals = [lambda * (x^2)^0.926 for x in q_points]
    df_evals = [lambda * 1.852*x / (x^2)^0.074 for x in q_points]
    df_evals[isnan.(df_evals)] = 0.0 # This will only affect q = 0.
    return [f_evals[i] + (q - q_points[i]) * df_evals[i] for i in 1:n]
end

"Return values that approximate the Darcy-Weisbach head loss constraint."
function construct_dw_separators(q::JuMP.Variable, lambda::Float64, n::Int = 3)
    q_points = linspace(getlowerbound(q), getupperbound(q), n)
    f_evals = [lambda * x^2 for x in q_points]
    df_evals = [2.0 * lambda * x for x in q_points]
    return [f_evals[i] + (q - q_points[i]) * df_evals[i] for i in 1:n]
end

"Piecewise linear Darcy-Weisbach constraint with unknown direction variables."
function constraint_dw_unknown_direction{T <: StandardMILPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Compute the relevant piecewise linear data.
    breakpoints = linspace(getlowerbound(q), getupperbound(q), 50)
    rhs = piecewiselinear(wm.model, q, breakpoints, (u) -> lambda * u * abs(u))

    # Add the piecewise linear constraint.
    @constraint(wm.model, h_i - h_j == rhs)
end

"Piecewise linear Darcy-Weisbach constraint with known direction variables."
function constraint_dw_known_direction{T <: StandardMILPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Compute the relevant piecewise linear data.
    breakpoints = linspace(getlowerbound(q), getupperbound(q), 50)
    rhs = piecewiselinear(wm.model, q, breakpoints, (u) -> lambda * u * abs(u))

    # Add the piecewise linear constraint.
    @constraint(wm.model, h_i - h_j == rhs)
end

"Relaxed piecewise linear Darcy-Weisbach constraint with unknown direction variables."
function constraint_dw_unknown_direction{T <: RelaxedMILPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)
    gamma = wm.var[:nw][n][:gamma][a]

    # Add constraints required to define gamma.
    constraint_define_gamma(wm, a, n)

    # Use the piecewise linear outer approximation.
    for cut in construct_dw_separators(q, lambda)
        @constraint(wm.model, gamma >= cut)
    end
end

"Relaxed piecewise linear Darcy-Weisbach constraint with known direction variables."
function constraint_dw_known_direction{T <: RelaxedMILPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Use the piecewise linear outer approximation.
    for cut in construct_dw_separators(q, lambda)
        @constraint(wm.model, dir * (h_i - h_j) >= cut)
    end
end

"Piecewise linear Hazen-Williams constraint with unknown direction variables."
function constraint_hw_unknown_direction{T <: StandardMILPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Compute the relevant piecewise linear data.
    breakpoints = linspace(getlowerbound(q), getupperbound(q), 50)
    rhs = piecewiselinear(wm.model, q, breakpoints, (u) -> lambda * u * (u^2)^0.426)

    # Add the piecewise linear constraint.
    @constraint(wm.model, h_i - h_j == rhs)
end

"Piecewise linear Hazen-Williams constraint with known direction variables."
function constraint_hw_known_direction{T <: StandardMILPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Compute the relevant piecewise linear data.
    breakpoints = linspace(getlowerbound(q), getupperbound(q), 50)
    rhs = piecewiselinear(wm.model, q, breakpoints, (u) -> lambda * u * (u^2)^0.426)

    # Add the piecewise linear constraint.
    @constraint(wm.model, h_i - h_j == rhs)
end

"Relaxed piecewise linear Hazen-Williams constraint with unknown direction variables."
function constraint_hw_unknown_direction{T <: RelaxedMILPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)
    gamma = wm.var[:nw][n][:gamma][a]

    # Add constraints required to define gamma.
    constraint_define_gamma(wm, a, n)

    # Use the piecewise linear outer approximation.
    for cut in construct_hw_separators(q, lambda)
        @constraint(wm.model, gamma >= cut)
    end
end

"Relaxed piecewise linear Hazen-Williams constraint with known direction variables."
function constraint_hw_known_direction{T <: RelaxedMILPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Use the piecewise linear outer approximation.
    for cut in construct_hw_separators(q, lambda)
        @constraint(wm.model, dir * (h_i - h_j) >= cut)
    end
end
