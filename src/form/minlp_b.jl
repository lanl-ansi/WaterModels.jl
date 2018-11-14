# Define MINLP-B implementations of water distribution models.

export MINLPBWaterModel, StandardMINLPBForm

abstract type AbstractMINLPBForm <: AbstractWaterFormulation end
abstract type StandardMINLPBForm <: AbstractMINLPBForm end

"The default MINLP-B model."
const MINLPBWaterModel = GenericWaterModel{StandardMINLPBForm}

"Default MINLP-B constructor."
MINLPBWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMINLPBForm; kwargs...)

"Non-convex Darcy-Weisbach constraint with unknown direction variables."
function constraint_dw_unknown_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: AbstractMINLPBForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Add a nonlinear constraint for the head loss.
    @NLconstraint(wm.model, h_i - h_j == lambda * q * abs(q))
end

"Quintic approximation of Bragalli, D'Ambrosio, Lee, Lodi, and Toth.
(http://www.optimization-online.org/DB_FILE/2008/03/1924.pdf)"
function hw_quintic_approx(model::JuMP.Model, q::JuMP.Variable, delta::Float64 = 0.1)
    c_5 = 3.0 * delta^(1.852 - 5) / 8.0 +
          1.0 / 8.0 * (1.852 - 1) * 1.852 * delta^(1.852 - 5) -
          3.0 / 8.0 * 1.852 * delta^(1.852 - 5)

    c_3 = -5.0 * delta^(1.852 - 3) / 4.0 -
          1.0 / 4.0 * (1.852 - 1) * 1.852 * delta^(1.852 - 3) +
          5.0 / 4.0 * 1.852 * delta^(1.852 - 3)

    c_1 = 15.0 * delta^(1.852 - 1) / 8.0 +
          1.0 / 8.0 * (1.852 - 1) * 1.852 * delta^(1.852 - 1) -
          7.0 / 8.0 * 1.852 * delta^(1.852 - 1)

    return @NLexpression(model, c_5 * q^5 + c_3 * q^3 + c_1 * q)
end

"Non-convex Hazen-Williams constraint for flow with unknown direction."
function constraint_hw_unknown_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: AbstractMINLPBForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Add big-M constraints for the piecewise indicator variable z.
    z = @variable(wm.model, category = :Bin)
    delta = @variable(wm.model, category = :Bin)
    M = 0.1 + max(abs(getlowerbound(q)), abs(getupperbound(q)))
    @constraint(wm.model, q <= -0.1 + M * delta + M * z)
    @constraint(wm.model, q >= 0.1 - M * (1 - delta) - M * z)

    # Add a piecewise nonlinear constraint for the head loss.
    inner_hw = hw_quintic_approx(wm.model, q, 0.1)
    inner_piece = @NLexpression(wm.model, z * lambda * inner_hw)

    # Use this for the time being.
    outer_piece = @NLexpression(wm.model, (1 - z) * lambda * q * (q^2)^0.426)
    @NLconstraint(wm.model, h_i - h_j == inner_piece + outer_piece)
end

"Non-convex Hazen-Williams constraint for flow with unknown direction."
function constraint_hw_unknown_direction_ne(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: AbstractMINLPBForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

    # Add big-M constraints for the piecewise indicator variable z.
    z = @variable(wm.model, category = :Bin)
    delta = @variable(wm.model, category = :Bin)
    M = 0.1 + max(abs(getlowerbound(q)), abs(getupperbound(q)))
    @constraint(wm.model, q <= -0.1 + M * delta + M * z)
    @constraint(wm.model, q >= 0.1 - M * (1 - delta) - M * z)

    # Add a piecewise nonlinear constraint for the head loss.
    inner_hw = hw_quintic_approx(wm.model, q, 0.1)
    inner_piece = @NLexpression(wm.model, z * inner_hw)

    # Add constraints required to define gamma.
    constraint_define_gamma_hw_ne(wm, a, n)

    # Define an auxiliary variable for the sum of the gamma variables.
    gamma_sum = @variable(wm.model, basename = "gamma_sum_$(n)_$(a)", start = 0)
    @constraint(wm.model, gamma_sum == sum(wm.var[:nw][n][:gamma][a]))

    # Use this for the time being.
    outer_piece = @NLexpression(wm.model, (1 - z) * q * (q^2)^0.426)
    @NLconstraint(wm.model, gamma_sum == inner_piece + outer_piece)
end
