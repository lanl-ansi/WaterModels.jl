# Defines MICP implementations of water distribution models.
export MICPWaterModel, StandardMICPForm

abstract type AbstractMICPForm <: AbstractWaterFormulation end
abstract type StandardMICPForm <: AbstractMICPForm end

"Default MICP (mixed-integer convex programming) model."
const MICPWaterModel = GenericWaterModel{StandardMICPForm}

"Default MICP constructor."
MICPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMICPForm; kwargs...)

function constraint_potential_loss(wm::GenericWaterModel{T}, a::Int, n_n::Int = wm.cnw) where T <: StandardMICPForm
    if !haskey(wm.con[:nw][n_n], :potential_loss_1)
        function_head_loss_hw(wm) # Register the head loss JuMP function.
        wm.con[:nw][n_n][:potential_loss_1] = Dict{Int, Dict{Int, ConstraintRef}}()
        wm.con[:nw][n_n][:potential_loss_2] = Dict{Int, Dict{Int, ConstraintRef}}()
    end

    wm.con[:nw][n_n][:potential_loss_1][a] = Dict{Int, ConstraintRef}()
    wm.con[:nw][n_n][:potential_loss_2][a] = Dict{Int, ConstraintRef}()

    dhp = wm.var[:nw][n_n][:dhp][a]
    dhn = wm.var[:nw][n_n][:dhn][a]
    L = wm.ref[:nw][n_n][:connection][a]["length"]

    for r in 1:length(wm.ref[:nw][n_n][:resistance][a])
        q_n_r = wm.var[:nw][n_n][:qn][a][r]
        q_p_r = wm.var[:nw][n_n][:qp][a][r]
        resistance = wm.ref[:nw][n_n][:resistance][a][r]

        con_1 = @NLconstraint(wm.model, dhp / L >= resistance * head_loss_hw(q_p_r))
        wm.con[:nw][n_n][:potential_loss_1][a][r] = con_1

        con_2 = @NLconstraint(wm.model, dhn / L >= resistance * head_loss_hw(q_n_r))
        wm.con[:nw][n_n][:potential_loss_2][a][r] = con_2
    end
end
