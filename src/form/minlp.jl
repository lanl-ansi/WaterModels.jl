# Defines MINLP implementations of water distribution models.
export MINLPWaterModel, StandardMINLPForm

abstract type AbstractMINLPForm <: AbstractWaterFormulation end
abstract type StandardMINLPForm <: AbstractMINLPForm end

"Default (convex) MINLP model."
const MINLPWaterModel = GenericWaterModel{StandardMINLPForm}

"Default MINLP constructor."
MINLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMINLPForm; kwargs...)

function constraint_potential_loss(wm::GenericWaterModel{T}, a::Int, n_n::Int = wm.cnw) where T <: StandardMINLPForm
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
