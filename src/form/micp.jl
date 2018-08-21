# Define MICP (mixed-integer convex program) implementations of water distribution models.

export MICPWaterModel, StandardMICPForm

""
@compat abstract type AbstractMICPForm <: AbstractMINLPForm end

""
@compat abstract type StandardMICPForm <: AbstractMICPForm end

"The default MICP (mixed-integer convex program) model is a relaxation of the non-convex MINLP model."
const MICPWaterModel = GenericWaterModel{StandardMICPForm}

"Default MICP constructor."
MICPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMICPForm; kwargs...)

"Create variables associated with the head for the MICP problem."
function variable_head{T <: AbstractMICPForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
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

"Constraints used to define the head difference in the MICP."
function constraint_define_gamma{T <: AbstractMICPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
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

"Convex (relaxed) Darcy-Weisbach constraint with known direction variables."
function constraint_dw_known_direction{T <: AbstractMICPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a convex quadratic constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) >= lambda * q^2)
end

"Convex (relaxed) Darcy-Weisbach constraint with unknown direction variables."
function constraint_dw_unknown_direction{T <: AbstractMICPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)
    gamma = wm.var[:nw][n][:gamma][a]

    # Add constraints required to define gamma.
    constraint_define_gamma(wm, a, n)

    # Add a convex quadratic constraint for the head loss.
    @NLconstraint(wm.model, gamma >= lambda * q^2)
end

"Convex (relaxed) Hazen-Williams constraint for flow with known direction."
function constraint_hw_known_direction{T <: AbstractMICPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a nonlinear constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) >= lambda * (q^2 + 1.0e-4)^0.926)
end

"Convex (relaxed) Hazen-Williams constraint for flow with unknown direction."
function constraint_hw_unknown_direction{T <: AbstractMICPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)
    gamma = wm.var[:nw][n][:gamma][a]

    # Add constraints required to define gamma.
    constraint_define_gamma(wm, a, n)

    # Add a nonlinear constraint for the head loss.
    @NLconstraint(wm.model, gamma >= lambda * (q^2 + 1.0e-4)^0.926)
end
