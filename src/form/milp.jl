# Define MILP implementations of water distribution models.

using PiecewiseLinearOpt

export MILPWaterModel, StandardMILPForm

"AbstractMILPForm is derived from AbstractMINLPForm"
@compat abstract type AbstractMILPForm <: AbstractWaterFormulation end

"StandardMILPForm is derived from StandardMINLPForm"
@compat abstract type StandardMILPForm <: AbstractMILPForm end

"The default MILP model is assumed to meet all constraints with equality."
const MILPWaterModel = GenericWaterModel{StandardMILPForm}

"Default MILP constructor (assumes equality constraints)."
MILPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMILPForm; kwargs...)

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
