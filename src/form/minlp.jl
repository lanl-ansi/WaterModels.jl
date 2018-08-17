# Define MINLP implementations of water distribution models.

export MINLPWaterModel, StandardMINLPForm

""
@compat abstract type AbstractMINLPForm <: AbstractWaterFormulation end

""
@compat abstract type StandardMINLPForm <: AbstractMINLPForm end

const MINLPWaterModel = GenericWaterModel{StandardMINLPForm}

"Default MINLP constructor."
MINLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMINLPForm; kwargs...)

"Darcy-Weisbach constraint with unknown direction variables."
function constraint_dw_unknown_direction{T <: AbstractMINLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    gamma = wm.var[:nw][n][:gamma][a]
    q = wm.var[:nw][n][:q][a]
    viscosity = wm.ref[:nw][n][:options]["viscosity"]
    lambda = calc_friction_factor_dw(wm.ref[:nw][n][:pipes][a], viscosity)
    @constraint(wm.model, gamma == lambda * q^2)
end

# From http://www.optimization-online.org/DB_FILE/2008/03/1924.pdf
function hw_approx(model, q, delta::Float64 = 0.1)
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

"Hazen-Williams constraint for flow with unknown direction."
function constraint_hw_unknown_direction{T <: AbstractMINLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Get source and target nodes corresponding to the edge.
    i = wm.ref[:nw][n][:pipes][a]["node1"]
    j = wm.ref[:nw][n][:pipes][a]["node2"]

    # Collect variables needed for the constraint.
    q = wm.var[:nw][n][:q][a]
    h_i = wm.var[:nw][n][:h][i]
    h_j = wm.var[:nw][n][:h][j]
    lambda = calc_friction_factor_hw(wm.ref[:nw][n][:pipes][a])

    # Add constraints for the piecewise indicator variable z.
    z = @variable(wm.model, category = :Bin)
    delta = @variable(wm.model, category = :Bin)
    M = max(abs(getupperbound(q)), abs(getlowerbound(q)))
    @constraint(wm.model, q <= -1.0e-3 + M * delta + M * z)
    @constraint(wm.model, q >= 1.0e-3 - M * (1 - delta) - M * z)

    # Add a piecewise nonlinear constraint for the head loss.
    inner_hw = hw_approx(wm.model, q, 1.0e-3)
    inner_piece = @NLexpression(wm.model, lambda * inner_hw)
    outer_piece = @NLexpression(wm.model, lambda * q * (q^2 + 1.0e-7)^0.426)
    @NLconstraint(wm.model, h_i - h_j == z * inner_piece + (1 - z) * outer_piece)
end

"Hazen-Williams constraint for flow with known direction."
function constraint_hw_known_direction{T <: AbstractMINLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Get source and target nodes corresponding to the edge.
    i = wm.ref[:nw][n][:pipes][a]["node1"]
    j = wm.ref[:nw][n][:pipes][a]["node2"]

    # Get the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])

    # Collect variables needed for the constraint.
    q = wm.var[:nw][n][:q][a]
    h_i = wm.var[:nw][n][:h][i]
    h_j = wm.var[:nw][n][:h][j]
    lambda = calc_friction_factor_hw(wm.ref[:nw][n][:pipes][a])

    # Fix the direction of the flow.
    setlowerbound(q, dir == 1 ? 0.0 : getlowerbound(q))
    setupperbound(q, dir == 1 ? getupperbound(q) : 0.0)

    # Add a nonlinear constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) == (q^2 + 1.0e-7)^0.926)
end
