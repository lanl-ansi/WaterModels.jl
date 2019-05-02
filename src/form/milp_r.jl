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
function construct_hw_separators(q::JuMP.VariableRef, lambda::Float64, n::Int = 3)
    q_points = LinRange(lower_bound(q), upper_bound(q), n)
    f_evals = [lambda * (x^2)^0.926 for x in q_points]
    df_evals = [lambda * 1.852*x / (x^2)^0.074 for x in q_points]
    nan_indices = findall(isnan, df_evals)
    setindex!(df_evals, zeros(size(nan_indices, 1)), nan_indices)
    return [f_evals[i] + (q - q_points[i]) * df_evals[i] for i in 1:n]
end

"Return values that approximate the Darcy-Weisbach head loss constraint."
function construct_dw_separators(q::JuMP.VariableRef, lambda::Float64, n::Int = 3)
    q_points = LinRange(lower_bound(q), upper_bound(q), n)
    f_evals = [lambda * x^2 for x in q_points]
    df_evals = [2.0 * lambda * x for x in q_points]
    return [f_evals[i] + (q - q_points[i]) * df_evals[i] for i in 1:n]
end

function constraint_potential_loss(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: StandardMILPRForm
    if !haskey(wm.con[:nw][n], :potential_loss_1)
        wm.con[:nw][n][:potential_loss_1] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        wm.con[:nw][n][:potential_loss_2] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    wm.con[:nw][n][:potential_loss_1][a] = Dict{Int, JuMP.ConstraintRef}()
    wm.con[:nw][n][:potential_loss_2][a] = Dict{Int, JuMP.ConstraintRef}()

    dhp = wm.var[:nw][n][:dhp][a]
    dhn = wm.var[:nw][n][:dhn][a]
    L = wm.ref[:nw][n][:links][a]["length"]

    # TODO: Not efficient... we need another method for storing resistances.
    R = calc_resistances_hw(wm, n)

    #for r in keys(wm.var[:nw][n][:xr][a])
    #    q_n_r = wm.var[:nw][n][:qn][a][r]
    #    q_p_r = wm.var[:nw][n][:qp][a][r]

    #    con_1 = @NLconstraint(wm.model, dhp / L >= R[a][r] * head_loss_hw_sep(q_p_r))
    #    wm.con[:nw][n][:potential_loss_1][a][r] = con_1

    #    con_2 = @NLconstraint(wm.model, dhn / L >= R[a][r] * head_loss_hw_sep(q_n_r))
    #    wm.con[:nw][n][:potential_loss_2][a][r] = con_2
    #end
end
