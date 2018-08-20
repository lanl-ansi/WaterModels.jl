# Define non-convex MINLP implementations of water distribution models.

export MINLPWaterModel, ConvexMINLPWaterModel, StandardMINLPForm, ConvexMINLPForm

""
@compat abstract type AbstractMINLPForm <: AbstractWaterFormulation end

""
@compat abstract type StandardMINLPForm <: AbstractMINLPForm end

""
@compat abstract type ConvexMINLPForm <: AbstractMINLPForm end

"The default MINLP model is assumed to be non-convex."
const MINLPWaterModel = GenericWaterModel{StandardMINLPForm}

"The convex (relaxed) MINLP has an explicit 'Convex' label."
const ConvexMINLPWaterModel = GenericWaterModel{ConvexMINLPForm}

"Default MINLP constructor (assumes the non-convex form)."
MINLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMINLPForm; kwargs...)

"Convex MINLP constructor."
ConvexMINLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, ConvexMINLPForm; kwargs...)

"Set new bounds for q given some specified direction of flow (-1 or 1)."
function fix_flow_direction(q::JuMP.Variable, direction::Int)
    # Fix the direction of the flow.
    setlowerbound(q, direction == 1 ? 0.0 : getlowerbound(q))
    setupperbound(q, direction == 1 ? getupperbound(q) : 0.0)
end

"Get variables commonly used in the construction of head loss constraints."
function get_common_variables{T <: AbstractMINLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Get source and target nodes corresponding to the edge.
    i = wm.ref[:nw][n][:pipes][a]["node1"]
    j = wm.ref[:nw][n][:pipes][a]["node2"]

    # Collect variables needed for the constraint.
    q = wm.var[:nw][n][:q][a]
    h_i = wm.var[:nw][n][:h][i]
    h_j = wm.var[:nw][n][:h][j]

    # Return the variables.
    return q, h_i, h_j
end

"Create variables associated with the head for the convex MINLP problem."
function variable_head{T <: ConvexMINLPForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    # Create head variables.
    variable_head_common(wm, n)

    # Create variables that correspond to the absolute value of the head difference.
    diff_min, diff_max = calc_head_difference_bounds(wm.ref[:nw][n][:pipes])
    wm.var[:nw][n][:gamma] = @variable(wm.model, [id in keys(wm.ref[:nw][n][:pipes])],
                                       lowerbound = diff_min[id], upperbound = diff_max[id],
                                       basename = "gamma_$(n)")

    # Create variables that correspond to flow moving from i to j.
    wm.var[:nw][n][:yp] = @variable(wm.model, [id in keys(wm.ref[:nw][n][:pipes])],
                                    category = :Bin, basename = "yp_$(n)")

    # Create variables that correspond to flow moving from j to i.
    wm.var[:nw][n][:yn] = @variable(wm.model, [id in keys(wm.ref[:nw][n][:pipes])],
                                    category = :Bin, basename = "yn_$(n)")
end

"Create variables associated with the head for the nonconvex MINLP problem."
function variable_head{T <: StandardMINLPForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    variable_head_common(wm, n)
end

"Get variables and constants used in the construction of Darcy-Weisbach constraints."
function get_dw_requirements{T <: AbstractMINLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    q, h_i, h_j = get_common_variables(wm, a, n)
    viscosity = wm.ref[:nw][n][:options]["viscosity"]
    lambda = calc_friction_factor_dw(wm.ref[:nw][n][:pipes][a], viscosity)
    return q, h_i, h_j, viscosity, lambda
end

"Non-convex Darcy-Weisbach constraint with known direction variables."
function constraint_dw_known_direction{T <: StandardMINLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a non-convex quadratic constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) == lambda * q^2)
end

"Convex (relaxed) Darcy-Weisbach constraint with known direction variables."
function constraint_dw_known_direction{T <: ConvexMINLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a convex quadratic constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) >= lambda * q^2)
end

"Non-convex Darcy-Weisbach constraint with unknown direction variables."
function constraint_dw_unknown_direction{T <: StandardMINLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Add a nonlinear constraint for the head loss.
    @NLconstraint(wm.model, h_i - h_j == lambda * q * abs(q))
end

"Convex (relaxed) Darcy-Weisbach constraint with unknown direction variables."
function constraint_dw_unknown_direction{T <: ConvexMINLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)
    gamma = wm.var[:nw][n][:gamma][a]

    # Add constraints required to define gamma.
    constraint_define_gamma(wm, a, n)

    # Add a convex quadratic constraint for the head loss.
    @NLconstraint(wm.model, gamma >= lambda * q^2)
end

"Get variables and constants used in the construction of Hazen-Williams constraints."
function get_hw_requirements{T <: AbstractMINLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    q, h_i, h_j = get_common_variables(wm, a, n)
    lambda = calc_friction_factor_hw(wm.ref[:nw][n][:pipes][a])
    return q, h_i, h_j, lambda
end

"Non-convex Hazen-Williams constraint for flow with known direction."
function constraint_hw_known_direction{T <: StandardMINLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a non-convex constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) == (q^2 + 1.0e-4)^0.926)
end

"Convex (relaxed) Hazen-Williams constraint for flow with known direction."
function constraint_hw_known_direction{T <: ConvexMINLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a nonlinear constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) >= (q^2 + 1.0e-4)^0.926)
end

"Quintic approximation from http://www.optimization-online.org/DB_FILE/2008/03/1924.pdf"
function hw_approx(model::JuMP.Model, q::JuMP.Variable, delta::Float64 = 0.1)
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
function constraint_hw_unknown_direction{T <: StandardMINLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Add big-M constraints for the piecewise indicator variable z.
    z = @variable(wm.model, category = :Bin)
    delta = @variable(wm.model, category = :Bin)
    M = 0.1 + max(abs(getlowerbound(q)), abs(getupperbound(q)))
    @constraint(wm.model, q <= -0.1 + M * delta + M * z)
    @constraint(wm.model, q >= 0.1 - M * (1 - delta) - M * z)

    # Add a piecewise nonlinear constraint for the head loss.
    inner_hw = hw_approx(wm.model, q, 0.1)
    inner_piece = @NLexpression(wm.model, z * lambda * inner_hw)

    # TODO: Use this when we finally have a working non-convex solver interface
    # that allows for user-defined functions and their derivatives.
    # outer_piece = @NLexpression(wm.model, (1 - z) * lambda * hw_q(q))

    # Use this for the time being.
    outer_piece = @NLexpression(wm.model, (1 - z) * lambda * q * (q^2 + 1.0e-4)^0.426)
    @NLconstraint(wm.model, h_i - h_j == inner_piece + outer_piece)
end

"Constraints used to define the head difference in the convex (relaxed) MINLP."
function constraint_define_gamma{T <: ConvexMINLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Get source and target nodes corresponding to the edge.
    i = wm.ref[:nw][n][:pipes][a]["node1"]
    j = wm.ref[:nw][n][:pipes][a]["node2"]

    # Collect variables needed for the constraint.
    q = wm.var[:nw][n][:q][a]
    h_i = wm.var[:nw][n][:h][i]
    h_j = wm.var[:nw][n][:h][j]

    # Collect other variables needed for the constraint.
    gamma = wm.var[:nw][n][:gamma][a]
    y_p = wm.var[:nw][n][:yp][a]
    y_n = wm.var[:nw][n][:yn][a]

    # Get additional data related to the variables.
    h_i_lb = getlowerbound(h_i)
    h_i_ub = getupperbound(h_i)
    h_j_lb = getlowerbound(h_j)
    h_j_ub = getupperbound(h_j)

    # Add the required constraints to define gamma.
    @constraint(wm.model, h_j - h_i + (h_i_lb - h_j_ub) * (y_p - y_n + 1) <= gamma)
    @constraint(wm.model, h_i - h_j + (h_i_ub - h_j_lb) * (y_p - y_n - 1) <= gamma)
    @constraint(wm.model, h_j - h_i + (h_i_ub - h_j_lb) * (y_p - y_n + 1) >= gamma)
    @constraint(wm.model, h_i - h_j + (h_i_lb - h_j_ub) * (y_p - y_n - 1) >= gamma)

    # Get the sum of all junction demands.
    junctions = values(wm.ref[:nw][n][:junctions])
    sum_demand = sum(junction["demand"] for junction in junctions)

    # Add the required constraints to define y_p and y_n.
    @constraint(wm.model, (y_p - 1) * sum_demand <= q)
    @constraint(wm.model, (1 - y_n) * sum_demand >= q)
    @constraint(wm.model, (1 - y_p) * (h_i_lb - h_j_ub) <= h_i - h_j)
    @constraint(wm.model, (1 - y_n) * (h_i_ub - h_j_lb) >= h_i - h_j)
    @constraint(wm.model, y_p + y_n == 1)
end

"Convex (relaxed) Hazen-Williams constraint for flow with unknown direction."
function constraint_hw_unknown_direction{T <: ConvexMINLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)
    gamma = wm.var[:nw][n][:gamma][a]

    # Add constraints required to define gamma.
    constraint_define_gamma(wm, a, n)

    # Add a nonlinear constraint for the head loss.
    @NLconstraint(wm.model, gamma >= (q^2 + 1.0e-4)^0.926)
end
