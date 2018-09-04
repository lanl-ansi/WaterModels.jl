"Set new bounds for q given some specified direction of flow (-1 or 1)."
function fix_flow_direction(q::JuMP.Variable, direction::Int)
    # Fix the direction of the flow.
    setlowerbound(q, direction == 1 ? 0.0 : getlowerbound(q))
    setupperbound(q, direction == 1 ? getupperbound(q) : 0.0)
end

"Get variables commonly used in the construction of head loss constraints."
function get_common_variables{T <: AbstractWaterFormulation}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
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

"Get variables and constants used in the construction of Darcy-Weisbach constraints."
function get_dw_requirements{T <: AbstractWaterFormulation}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    q, h_i, h_j = get_common_variables(wm, a, n)
    viscosity = wm.ref[:nw][n][:options]["viscosity"]
    lambda = calc_friction_factor_dw(wm.ref[:nw][n][:pipes][a], viscosity)
    return q, h_i, h_j, viscosity, lambda
end

"Get variables and constants used in the construction of Hazen-Williams constraints."
function get_hw_requirements{T <: AbstractWaterFormulation}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    q, h_i, h_j = get_common_variables(wm, a, n)

    lambda = nothing
    if haskey(wm.var[:nw][n][:lambda], a)
        lambda = wm.var[:nw][n][:lambda][a]
    else
        lambda = calc_friction_factor_hw(wm.ref[:nw][n][:pipes][a])
    end

    return q, h_i, h_j, lambda
end


"These problem forms use binary variables to specify flow direction."
AbstractRelaxedForm = Union{AbstractMICPForm, AbstractMILPRForm}

"Create variables associated with the head for the MICP and MILP-R problems."
function variable_head{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
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

"Constraints used to define the head difference in the MICP and MILP-R problems."
function constraint_define_gamma{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

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

"These problem forms use exact versions of the head loss equation when the flow direction is fixed."
AbstractUndirectedForm = Union{AbstractMILPForm, AbstractMINLPBForm, AbstractNLPForm}

"Create variables associated with the head for forms of the problem without direction variables."
function variable_head{T <: AbstractUndirectedForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    variable_head_common(wm, n)
end

"These problem forms use exact versions of the head loss equation when the flow direction is fixed."
AbstractExactForm = Union{AbstractMINLPBForm, AbstractNLPForm}

"Exact Hazen-Williams constraint for flow with known direction."
function constraint_hw_known_direction{T <: AbstractExactForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a non-convex constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) == lambda * (q^2)^0.926)
end

"Exact Darcy-Weisbach constraint with known direction."
function constraint_dw_known_direction{T <: AbstractExactForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a non-convex quadratic constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) == lambda * q^2)
end
