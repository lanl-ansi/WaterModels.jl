# Define MICP (mixed-integer convex program) implementations of water distribution models.
export MICPWaterModel, StandardMICPForm

abstract type AbstractMICPForm <: AbstractWaterFormulation end
abstract type StandardMICPForm <: AbstractMICPForm end

"The default MICP (mixed-integer convex program) model is a relaxation of the non-convex MINLP model."
const MICPWaterModel = GenericWaterModel{StandardMICPForm}

"Default MICP constructor."
MICPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMICPForm; kwargs...)

function variable_head(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852) where T <: AbstractMICPForm
    variable_pressure_head(wm, n)
    variable_directed_head_difference(wm, n)
end

function variable_flow(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852) where T <: AbstractMICPForm
    variable_undirected_flow(wm, n, bounded=true)
    variable_directed_flow(wm, n, alpha=alpha, bounded=true)
    variable_flow_direction(wm, n)
end

function variable_flow_ne(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852) where T <: AbstractMICPForm
    variable_directed_flow_ne(wm, n, alpha=alpha, bounded=true)
end

function constraint_resistance_selection_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractMICPForm
    constraint_directed_resistance_selection_ne(wm, a, n)
end

function constraint_potential_loss(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw; alpha::Float64=1.852) where T <: AbstractMICPForm
    constraint_directed_head_difference(wm, a, n)
    constraint_flow_direction_selection(wm, a, n)
    constraint_directed_potential_loss_ub(wm, a, n, alpha=alpha)
    constraint_directed_potential_loss(wm, a, n)
end

function constraint_potential_loss_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw; alpha::Float64=1.852) where T <: AbstractMICPForm
    constraint_flow_direction_selection_ne(wm, a, n)
    constraint_directed_potential_loss_ub_ne(wm, a, n, alpha=alpha)
    constraint_directed_potential_loss_ne(wm, a, n)
end

function constraint_flow_conservation(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractMICPForm
    constraint_directed_flow_conservation(wm, i, n)
end

function constraint_flow_conservation_ne(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractMICPForm
    constraint_directed_flow_conservation_ne(wm, i, n)
end

function constraint_link_flow_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractMICPForm
    constraint_link_directed_flow_ne(wm, a, n)
end

function constraint_link_flow(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractMICPForm
    constraint_link_directed_flow(wm, a, n)
end

function constraint_source_flow(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractMICPForm
    constraint_directed_source_flow(wm, i, n)
end

function constraint_sink_flow(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractMICPForm
    constraint_directed_sink_flow(wm, i, n)
end

function constraint_directed_potential_loss_ne(wm::GenericWaterModel{T}, a::Int, n::Int) where T <: AbstractMICPForm
    if !haskey(wm.con[:nw][n], :potential_loss_n_ne)
        wm.con[:nw][n][:potential_loss_n_ne] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        wm.con[:nw][n][:potential_loss_p_ne] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    wm.con[:nw][n][:potential_loss_n_ne][a] = Dict{Int, JuMP.ConstraintRef}()
    wm.con[:nw][n][:potential_loss_p_ne][a] = Dict{Int, JuMP.ConstraintRef}()

    L = wm.ref[:nw][n][:links][a]["length"]
    alpha = wm.ref[:nw][n][:alpha]

    for (r_id, r) in enumerate(wm.ref[:nw][n][:resistance][a])
        qn_ne = wm.var[:nw][n][:qn_ne][a][r_id]
        dhn = wm.var[:nw][n][:dhn][a]
        #con_n = JuMP.@NLconstraint(wm.model, r * f_alpha(qn_ne) - inv(L) * dhn <= 0.0)
        con_n = JuMP.@NLconstraint(wm.model, r * abs(qn_ne)^(1.852) - inv(L) * dhn <= 0.0)
        wm.con[:nw][n][:potential_loss_n_ne][a][r_id] = con_n

        qp_ne = wm.var[:nw][n][:qp_ne][a][r_id]
        dhp = wm.var[:nw][n][:dhp][a]
        #con_p = JuMP.@NLconstraint(wm.model, r * f_alpha(qp_ne) - inv(L) * dhp <= 0.0)
        con_p = JuMP.@NLconstraint(wm.model, r * abs(qp_ne)^(1.852) - inv(L) * dhp <= 0.0)
        wm.con[:nw][n][:potential_loss_p_ne][a][r_id] = con_p
    end
end

function constraint_directed_potential_loss(wm::GenericWaterModel{T}, a::Int, n::Int) where T <: AbstractMICPForm
    #if !haskey(wm.con[:nw][n], :potential_loss_n)
    #    wm.con[:nw][n][:potential_loss_n] = Dict{Int, JuMP.ConstraintRef}()
    #    wm.con[:nw][n][:potential_loss_p] = Dict{Int, JuMP.ConstraintRef}()
    #end

    #L = wm.ref[:nw][n][:links][a]["length"]
    #r = minimum(wm.ref[:nw][n][:resistance][a])

    #qn = wm.var[:nw][n][:qn][a]
    #dhn = wm.var[:nw][n][:dhn][a]
    #con_n = JuMP.@NLconstraint(wm.model, r * f_alpha(qn) - inv(L) * dhn <= 0.0)
    #wm.con[:nw][n][:potential_loss_n][a] = con_n

    #qp = wm.var[:nw][n][:qp][a]
    #dhp = wm.var[:nw][n][:dhp][a]
    #con_p = JuMP.@NLconstraint(wm.model, r * f_alpha(qp) - inv(L) * dhp <= 0.0)
    #wm.con[:nw][n][:potential_loss_p][a] = con_p
end

function objective_wf(wm::GenericWaterModel{T}, n::Int = wm.cnw) where T <: StandardMICPForm
    JuMP.set_objective_sense(wm.model, MOI.FEASIBILITY_SENSE)
end
