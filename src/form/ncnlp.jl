# Define NCNLP (mixed-integer convex program) implementations of water distribution models.

export NCNLPWaterModel, StandardNCNLPForm

abstract type AbstractNCNLPForm <: AbstractWaterFormulation end
abstract type StandardNCNLPForm <: AbstractNCNLPForm end

"The default NCNLP (mixed-integer convex program) model is a relaxation of the non-convex MINLP model."
const NCNLPWaterModel = GenericWaterModel{StandardNCNLPForm}

"Default NCNLP constructor."
NCNLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardNCNLPForm; kwargs...)

function variable_head(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852) where T <: AbstractNCNLPForm
    variable_pressure_head(wm, n)
end

function variable_flow(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852) where T <: AbstractNCNLPForm
    variable_undirected_flow(wm, n, alpha=alpha, bounded=true)
end

function variable_flow_ne(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852) where T <: AbstractNCNLPForm
    variable_undirected_flow_ne(wm, n, alpha=alpha, bounded=true)
end

function constraint_resistance_selection_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractNCNLPForm
    constraint_undirected_resistance_selection_ne(wm, a, n)
end

function constraint_potential_loss(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw; alpha::Float64=1.852) where T <: AbstractNCNLPForm
    constraint_undirected_potential_loss(wm, a, n)
end

function constraint_potential_loss_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw; alpha::Float64=1.852) where T <: AbstractNCNLPForm
    constraint_undirected_potential_loss_ne(wm, a, n)
end

function constraint_flow_conservation(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractNCNLPForm
    constraint_undirected_flow_conservation(wm, i, n)
end

function constraint_flow_conservation_ne(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractNCNLPForm
    constraint_undirected_flow_conservation_ne(wm, i, n)
end

function constraint_link_flow_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractNCNLPForm
    constraint_link_undirected_flow_ne(wm, a, n)
end

function constraint_link_flow(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractNCNLPForm
end

function constraint_source_flow(wm::GenericWaterModel, i::Int, n::Int=wm.cnw) where T <: AbstractNCNLPForm
    constraint_undirected_source_flow(wm, i, n)
end

function constraint_sink_flow(wm::GenericWaterModel, i::Int, n::Int=wm.cnw) where T <: AbstractNCNLPForm
    constraint_undirected_sink_flow(wm, i, n)
end

function constraint_directed_potential_loss_ne(wm::GenericWaterModel{T}, a::Int, n::Int) where T <: AbstractNCNLPForm
    if !haskey(wm.con[:nw][n], :potential_lossⁿᵉ⁻)
        wm.con[:nw][n][:potential_lossⁿᵉ⁻] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        wm.con[:nw][n][:potential_lossⁿᵉ⁺] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    wm.con[:nw][n][:potential_lossⁿᵉ⁻][a] = Dict{Int, JuMP.ConstraintRef}()
    wm.con[:nw][n][:potential_lossⁿᵉ⁺][a] = Dict{Int, JuMP.ConstraintRef}()

    L = wm.ref[:nw][n][:links][a]["length"]

    for (r_id, r) in enumerate(wm.ref[:nw][n][:resistance][a])
        qⁿᵉ⁻ = wm.var[:nw][n][:qⁿᵉ⁻][a][r_id]
        Δh⁻ = wm.var[:nw][n][:Δh⁻][a]
        con⁻ = JuMP.@NLconstraint(wm.model, r * f_alpha(qⁿᵉ⁻) - inv(L) * Δh⁻ <= 0.0)
        wm.con[:nw][n][:potential_lossⁿᵉ⁻][a][r_id] = con⁻

        qⁿᵉ⁺ = wm.var[:nw][n][:qⁿᵉ⁺][a][r_id]
        Δh⁺ = wm.var[:nw][n][:Δh⁺][a]
        con⁺ = JuMP.@NLconstraint(wm.model, r * f_alpha(qⁿᵉ⁺) - inv(L) * Δh⁺ <= 0.0)
        wm.con[:nw][n][:potential_lossⁿᵉ⁺][a][r_id] = con⁺
    end
end

function constraint_undirected_potential_loss(wm::GenericWaterModel{T}, a::Int, n::Int) where T <: AbstractNCNLPForm
    if !haskey(wm.con[:nw][n], :potential_loss)
        wm.con[:nw][n][:potential_loss] = Dict{Int, JuMP.ConstraintRef}()
    end

    if i in collect(ids(wm, n, :reservoirs))
        hᵢ = wm.ref[:nw][n][:reservoirs][i]["head"]
    else
        hᵢ = wm.var[:nw][n][:h][i]
    end

    j = wm.ref[:nw][n][:links][a]["node2"]

    if j in collect(ids(wm, n, :reservoirs))
        hⱼ = wm.ref[:nw][n][:reservoirs][j]["head"]
    else
        hⱼ = wm.var[:nw][n][:h][j]
    end

    L = wm.ref[:nw][n][:links][a]["length"]
    r = minimum(wm.ref[:nw][n][:resistance][a])
    q = wm.var[:nw][n][:q][a]

    con = JuMP.@NLconstraint(wm.model, r * f_alpha(q) - inv(L) * (hᵢ - hⱼ) == 0.0)
    wm.con[:nw][n][:potential_loss][a] = con
end
