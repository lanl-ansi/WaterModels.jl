# Define MILP-R implementations of water distribution models.

export MILPRWaterModel, StandardMILPRForm

"AbstractMILPRForm is derived from AbstractWaterFormulation"
abstract type AbstractMILPRForm <: AbstractWaterFormulation end

"StandardMILPRForm is derived from AbstractMILPRForm"
abstract type StandardMILPRForm <: AbstractMILPRForm end

"The MILP-R model relaxes the head loss constraint to an inequality."
const MILPRWaterModel = GenericWaterModel{StandardMILPRForm}

"MILP-R constructor."
MILPRWaterModel(data::Dict{String, Any}; kwargs...) = GenericWaterModel(data, StandardMILPRForm; kwargs...)

"Return values that approximate the Hazen-Williams head loss constraint."
function construct_hw_separators(q::JuMP.Variable, lambda::Float64, n::Int = 3)
    q_points = LinRange(getlowerbound(q), getupperbound(q), n)
    f_evals = [lambda * (x^2)^0.926 for x in q_points]
    df_evals = [lambda * 1.852*x / (x^2)^0.074 for x in q_points]
    nan_indices = findall(isnan, df_evals)
    setindex!(df_evals, zeros(size(nan_indices, 1)), nan_indices)
    return [f_evals[i] + (q - q_points[i]) * df_evals[i] for i in 1:n]
end

"Return values that approximate the Darcy-Weisbach head loss constraint."
function construct_dw_separators(q::JuMP.Variable, lambda::Float64, n::Int = 3)
    q_points = LinRange(getlowerbound(q), getupperbound(q), n)
    f_evals = [lambda * x^2 for x in q_points]
    df_evals = [2.0 * lambda * x for x in q_points]
    return [f_evals[i] + (q - q_points[i]) * df_evals[i] for i in 1:n]
end

"Piecewise linear Darcy-Weisbach inequality constraints with unknown direction variables."
function constraint_dw_unknown_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMILPRForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Add constraints to restrict the direction.
    constraint_restrict_direction(wm, a, n)

    # Add constraints required to define gamma.
    constraint_define_gamma(wm, a, n)

    # Use the piecewise linear outer approximation.
    gamma = wm.var[:nw][n][:gamma][a]
    for cut in construct_dw_separators(q, lambda)
        @constraint(wm.model, gamma >= cut)
    end
end

"Piecewise linear Hazen-Williams inequality constraints with unknown direction variables."
function constraint_hw_unknown_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMILPRForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Add constraints to restrict the direction.
    constraint_restrict_direction(wm, a, n)

    # Add constraints required to define gamma.
    constraint_define_gamma(wm, a, n)

    # Use the piecewise linear outer approximation.
    gamma = wm.var[:nw][n][:gamma][a]
    for cut in construct_hw_separators(q, 1.0, 21)
        @constraint(wm.model, gamma >= cut)
    end
end

"Piecewise linear Darcy-Weisbach inequality constraints with unknown direction variables."
function constraint_dw_unknown_direction_ne(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMILPRForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

    # Add constraints to restrict the direction.
    constraint_restrict_direction(wm, a, n)

    # Add constraints required to define gamma.
    constraint_define_gamma_dw_ne(wm, a, n)

    # Use the piecewise linear outer approximation.
    for cut in construct_dw_separators(q, 1.0)
        @constraint(wm.model, sum(wm.var[:nw][n][:gamma][a]) >= cut)
    end
end

"Piecewise linear Hazen-Williams inequality constraints with unknown direction variables."
function constraint_hw_unknown_direction_ne(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: StandardMILPRForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

    # Add constraints to restrict the direction.
    constraint_restrict_direction(wm, a, n)

    # Add constraints required to define gamma.
    constraint_define_gamma_hw_ne(wm, a, n)

    # Define an auxiliary variable for the sum of the gamma variables.
    gamma_sum = wm.var[:nw][n][:gamma_sum][a]
    @constraint(wm.model, gamma_sum == sum(wm.var[:nw][n][:gamma][a]))

    # Compute the number of separators assuming intervals of 0.01.
    q_diff = getupperbound(q) - getlowerbound(q)
    #num_separators = max(3, Int(1 + ceil(q_diff / 0.01)))
    num_separators = 250

    # Use the piecewise linear outer approximation.
    for cut in construct_hw_separators(q, 1.0, num_separators)
        @constraint(wm.model, gamma_sum >= cut)
    end
end

"Piecewise linear Darcy-Weisbach inequality constraints with known direction variables."
function constraint_dw_known_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMILPRForm
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

"Piecewise linear Hazen-Williams inequality constraints with known direction variables."
function constraint_hw_known_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMILPRForm
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
