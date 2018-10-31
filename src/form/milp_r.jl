# Define MILP-R implementations of water distribution models.

export MILPRWaterModel, StandardMILPRForm

"AbstractMILPRForm is derived from AbstractWaterFormulation"
@compat abstract type AbstractMILPRForm <: AbstractWaterFormulation end

"StandardMILPRForm is derived from AbstractMILPRForm"
@compat abstract type StandardMILPRForm <: AbstractMILPRForm end

"The MILP-R model relaxes the head loss constraint to an inequality."
const MILPRWaterModel = GenericWaterModel{StandardMILPRForm}

"MILP-R constructor."
MILPRWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMILPRForm; kwargs...)

"Return values that approximate the Hazen-Williams head loss constraint."
function construct_hw_separators(q::JuMP.Variable, lambda::Float64, n::Int = 3)
    q_points = linspace(getlowerbound(q), getupperbound(q), n)
    f_evals = [lambda * (x^2)^0.926 for x in q_points]
    df_evals = [lambda * 1.852*x / (x^2)^0.074 for x in q_points]
    df_evals[isnan.(df_evals)] = 0.0 # This will only affect the value at x = 0.
    return [f_evals[i] + (q - q_points[i]) * df_evals[i] for i in 1:n]
end

"Return values that approximate the Darcy-Weisbach head loss constraint."
function construct_dw_separators(q::JuMP.Variable, lambda::Float64, n::Int = 3)
    q_points = linspace(getlowerbound(q), getupperbound(q), n)
    f_evals = [lambda * x^2 for x in q_points]
    df_evals = [2.0 * lambda * x for x in q_points]
    return [f_evals[i] + (q - q_points[i]) * df_evals[i] for i in 1:n]
end

"Piecewise linear Darcy-Weisbach inequality constraints with unknown direction variables."
function constraint_dw_unknown_direction{T <: StandardMILPRForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
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

"Piecewise linear Darcy-Weisbach inequality constraints with unknown direction variables."
function constraint_dw_unknown_direction_ne{T <: StandardMILPRForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

    # Add constraints required to define gamma.
    constraint_define_gamma_dw_ne(wm, a, n)

    # Use the piecewise linear outer approximation.
    for cut in construct_dw_separators(q, 1.0)
        @constraint(wm.model, sum(wm.var[:nw][n][:gamma][a]) >= cut)
    end
end

"Piecewise linear Darcy-Weisbach inequality constraints with known direction variables."
function constraint_dw_known_direction{T <: StandardMILPRForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
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

"Piecewise linear Hazen-Williams inequality constraints with unknown direction variables."
function constraint_hw_unknown_direction{T <: StandardMILPRForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)
    gamma = wm.var[:nw][n][:gamma][a]

    # Add constraints required to define gamma.
    constraint_define_gamma(wm, a, n)

    # Use the piecewise linear outer approximation.
    for cut in construct_hw_separators(q, 1.0, 21)
        @constraint(wm.model, gamma >= cut)
    end
end

"Piecewise linear Hazen-Williams inequality constraints with unknown direction variables."
function constraint_hw_unknown_direction_ne{T <: StandardMILPRForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

    # Add constraints required to define gamma.
    constraint_define_gamma_hw_ne(wm, a, n)

    # Define an auxiliary variable for the sum of the gamma variables.
    gamma_sum = wm.var[:nw][n][:gamma_sum][a]
    @constraint(wm.model, gamma_sum == sum(wm.var[:nw][n][:gamma][a]))

    ## Use the piecewise linear outer approximation.
    #gamma_ub = getupperbound(h_i) - getlowerbound(h_j)
    #gamma_lb = getlowerbound(h_i) - getupperbound(h_j)
    #num_points = (gamma_ub - gamma_lb) / 0.01
    for cut in construct_hw_separators(q, 1.0, 20)
        @constraint(wm.model, gamma_sum >= cut)
    end
end

"Piecewise linear Hazen-Williams inequality constraints with known direction variables."
function constraint_hw_known_direction{T <: StandardMILPRForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
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
