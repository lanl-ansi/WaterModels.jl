# Define MILP implementations of water distribution models.

using PiecewiseLinearOpt

export MILPWaterModel, StandardMILPForm

"AbstractMILPForm is derived from AbstractMINLPForm"
abstract type AbstractMILPForm <: AbstractWaterFormulation end

"StandardMILPForm is derived from StandardMINLPForm"
abstract type StandardMILPForm <: AbstractMILPForm end

"The default MILP model is assumed to meet all constraints with equality."
const MILPWaterModel = GenericWaterModel{StandardMILPForm}

"Default MILP constructor (assumes equality constraints)."
MILPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMILPForm; kwargs...)

"Piecewise linear Darcy-Weisbach constraint with unknown direction variables."
function constraint_dw_unknown_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMILPForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Compute the relevant piecewise linear data.
    breakpoints = linspace(getlowerbound(q), getupperbound(q), 50)
    rhs = piecewiselinear(wm.model, q, breakpoints, (u) -> lambda * u * abs(u), method = :ZigZag)

    # Add the piecewise linear constraint.
    @constraint(wm.model, h_i - h_j == rhs)
end

"Piecewise linear Darcy-Weisbach constraint for flow with unknown direction."
function constraint_dw_unknown_direction_ne(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMILPForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

    # Add constraints required to define gamma.
    constraint_define_gamma_dw_ne(wm, a, n)

    # Compute the relevant piecewise linear data.
    breakpoints = linspace(getlowerbound(q), getupperbound(q), 50)
    rhs = piecewiselinear(wm.model, q, breakpoints, (u) -> u * abs(u), method = :ZigZag)

    # Add the piecewise linear constraint.
    @constraint(wm.model, sum(wm.var[:nw][n][:gamma][a]) == rhs)
end

"Piecewise linear Darcy-Weisbach constraint with known direction variables."
function constraint_dw_known_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMILPForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Compute the relevant piecewise linear data.
    breakpoints = linspace(getlowerbound(q), getupperbound(q), 50)
    rhs = piecewiselinear(wm.model, q, breakpoints, (u) -> lambda * u * abs(u), method = :ZigZag)

    # Add the piecewise linear constraint.
    @constraint(wm.model, h_i - h_j == rhs)
end

"Piecewise linear Hazen-Williams constraint with unknown direction variables."
function constraint_hw_unknown_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMILPForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Compute the relevant piecewise linear data.
    breakpoints = linspace(getlowerbound(q), getupperbound(q), 50)
    rhs = piecewiselinear(wm.model, q, breakpoints, (u) -> lambda * u * (u^2)^0.426, method = :ZigZag)

    # Add the piecewise linear constraint.
    @constraint(wm.model, h_i - h_j == rhs)
end

"Piecewise linear Hazen-Williams constraint for flow with unknown direction."
function constraint_hw_unknown_direction_ne(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMILPForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

    # Add constraints required to define gamma.
    constraint_define_gamma_hw_ne(wm, a, n)

    # Compute the relevant piecewise linear data.
    breakpoints = linspace(getlowerbound(q), getupperbound(q), 20)
    rhs = piecewiselinear(wm.model, q, breakpoints, (u) -> u * (u^2)^0.426, method = :ZigZag)

    # Add the piecewise linear constraint.
    @constraint(wm.model, sum(wm.var[:nw][n][:gamma][a]) == rhs)
end

"Piecewise linear Hazen-Williams constraint with known direction variables."
function constraint_hw_known_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMILPForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Compute the relevant piecewise linear data.
    breakpoints = linspace(getlowerbound(q), getupperbound(q), 21)
    rhs = piecewiselinear(wm.model, q, breakpoints, (u) -> lambda * u * (u^2)^0.426, method = :ZigZag)

    # Add the piecewise linear constraint.
    @constraint(wm.model, h_i - h_j == rhs)
end

"Piecewise linear Hazen-Williams constraint for flow with unknown direction."
function constraint_hw_known_direction_ne(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMILPForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

    # Add constraints required to define gamma.
    constraint_define_gamma_hw_ne(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Fix the head difference directions.
    for gamma_d in wm.var[:nw][n][:gamma][a]
        fix_flow_direction(gamma, dir)
    end

    # Compute the relevant piecewise linear data.
    breakpoints = linspace(getlowerbound(q), getupperbound(q), 50)
    rhs = piecewiselinear(wm.model, q, breakpoints, (u) -> u * (u^2)^0.426, method = :ZigZag)

    # Add the piecewise linear constraint.
    @constraint(wm.model, sum(wm.var[:nw][n][:gamma][a]) == rhs)
end
