# Define MINLP-B implementations of water distribution models.

export MINLPBWaterModel, StandardMINLPBForm

""
@compat abstract type AbstractMINLPBForm <: AbstractWaterFormulation end

""
@compat abstract type StandardMINLPBForm <: AbstractMINLPBForm end

"The default MINLP-B model."
const MINLPBWaterModel = GenericWaterModel{StandardMINLPBForm}

"Default MINLP-B constructor."
MINLPBWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMINLPBForm; kwargs...)

"Non-convex Darcy-Weisbach constraint with unknown direction variables."
function constraint_dw_unknown_direction{T <: AbstractMINLPBForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
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
function constraint_hw_unknown_direction{T <: AbstractMINLPBForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
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

    # TODO: Use this when we finally have a working non-convex solver interface
    # that allows for user-defined functions and their derivatives.
    # outer_piece = @NLexpression(wm.model, (1 - z) * lambda * hw_q(q))

    # Use this for the time being.
    outer_piece = @NLexpression(wm.model, (1 - z) * lambda * q * (q^2)^0.426)
    @NLconstraint(wm.model, h_i - h_j == inner_piece + outer_piece)
end
