# Define MILP implementations of water distribution models.

export MILPWaterModel, StandardMILPForm

""
@compat abstract type AbstractMILPForm <: AbstractWaterFormulation end

""
@compat abstract type StandardMILPForm <: AbstractMILPForm end

const MILPWaterModel = GenericWaterModel{StandardMILPForm}

"default MILP constructor"
MILPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMILPForm; kwargs...)

function construct_hw_separators(q::JuMP.Variable, lambda::Float64, n::Int = 3)
    q_points = linspace(getlowerbound(q), getupperbound(q), n)
    f_evals = [lambda * (x^2)^0.926 for x in q_points]
    df_evals = [lambda * 1.852*x / (x^2)^0.074 for x in q_points]
    df_evals[isnan.(df_evals)] = 0.0
    return [f_evals[i] + (q - q_points[i]) * df_evals[i] for i in 1:n]
end

function construct_dw_separators(q::JuMP.Variable, lambda::Float64, n::Int = 3)
    q_points = linspace(getlowerbound(q), getupperbound(q), n)
    f_evals = [lambda * x^2 for x in q_points]
    df_evals = [2 * lambda * x for x in q_points]
    return [f_evals[i] + (q - q_points[i]) * df_evals[i] for i in 1:n]
end

"Darcy-Weisbach constraint with unknown direction variables."
function constraint_dw_unknown_direction{T}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    gamma = wm.var[:nw][n][:gamma][a]
    q = wm.var[:nw][n][:q][a]
    viscosity = wm.ref[:nw][n][:options]["viscosity"]
    lambda = calc_friction_factor_dw(wm.ref[:nw][n][:pipes][a], viscosity)

    # Use the piecewise linear outer approximation.
    for cut in construct_dw_separators(q, lambda)
        @constraint(wm.model, gamma >= cut)
    end
end

"Hazen-Williams constraint with unknown direction variables."
function constraint_hw_unknown_direction{T}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    gamma = wm.var[:nw][n][:gamma][a]
    q = wm.var[:nw][n][:q][a]
    lambda = calc_friction_factor_hw(wm.ref[:nw][n][:pipes][a])

    # Use the piecewise linear outer approximation.
    for cut in construct_hw_separators(q, lambda)
        @constraint(wm.model, gamma >= cut)
    end
end
