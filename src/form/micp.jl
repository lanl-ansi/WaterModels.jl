# Define MICP (mixed-integer convex program) implementations of water distribution models.

function variable_head(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractMICPForm
    variable_pressure_head(wm, n)
    variable_directed_head_difference(wm, n)
end

function variable_flow(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractMICPForm
    variable_undirected_flow(wm, n, bounded=true)
    variable_directed_flow(wm, n, bounded=true)
    variable_flow_direction(wm, n)
end

function variable_flow_ne(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractMICPForm
    variable_directed_flow_ne(wm, n, bounded=true)
end

function variable_pump(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractMICPForm
end

function constraint_resistance_selection_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractMICPForm
    constraint_directed_resistance_selection_ne(wm, a, n)
end

function constraint_potential_loss_pump(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractMICPForm
end

# TODO see if this can be removed
function constraint_potential_loss_pipe_ne(wm::GenericWaterModel{T}, a::Int; nw::Int=wm.cnw, kwargs...) where T <: AbstractMICPForm
    constraint_potential_loss_pipe_ne(wm, nw, a)
    constraint_flow_direction_selection_ne(wm, a; nw=nw, kwargs...)
    constraint_potential_loss_ub_pipe_ne(wm, a; nw=nw, kwargs...)
end

function constraint_source_flow(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractMICPForm
    constraint_directed_source_flow(wm, i, n)
end

function constraint_sink_flow(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractMICPForm
    constraint_directed_sink_flow(wm, i, n)
end

function constraint_potential_loss_pipe_ne(wm::GenericWaterModel{T}, n::Int, a::Int) where T <: AbstractMICPForm
    if !haskey(con(wm, n), :potential_loss_n_ne)
        con(wm, n)[:potential_loss_n_ne] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        con(wm, n)[:potential_loss_p_ne] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    con(wm, n, :potential_loss_n_ne)[a] = Dict{Int, JuMP.ConstraintRef}()
    con(wm, n, :potential_loss_p_ne)[a] = Dict{Int, JuMP.ConstraintRef}()

    L = ref(wm, n, :links, a)["length"]

    for (r_id, r) in enumerate(ref(wm, n, :resistance, a))
        qn_ne = var(wm, n, :qn_ne, a)[r_id]
        dhn = var(wm, n, :dhn, a)
        con_n = JuMP.@NLconstraint(wm.model, r * f_alpha(qn_ne) - inv(L) * dhn <= 0.0)
        con(wm, n, :potential_loss_n_ne, a)[r_id] = con_n

        qp_ne = var(wm, n, :qp_ne, a)[r_id]
        dhp = var(wm, n, :dhp, a)
        con_p = JuMP.@NLconstraint(wm.model, r * f_alpha(qp_ne) - inv(L) * dhp <= 0.0)
        con(wm, n, :potential_loss_p_ne, a)[r_id] = con_p
    end
end

function constraint_potential_loss_pipe(wm::GenericWaterModel{T}, n::Int, a::Int) where T <: AbstractMICPForm
    if !haskey(con(wm, n), :potential_loss_n)
        con(wm, n)[:potential_loss_n] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:potential_loss_p] = Dict{Int, JuMP.ConstraintRef}()
    end

    L = ref(wm, n, :links, a)["length"]
    r = minimum(ref(wm, n, :resistance, a))

    qn = var(wm, n, :qn, a)
    dhn = var(wm, n, :dhn, a)
    con_n = JuMP.@NLconstraint(wm.model, r * f_alpha(qn) - inv(L) * dhn <= 0.0)
    con(wm, n, :potential_loss_n)[a] = con_n

    qp = var(wm, n, :qp, a)
    dhp = var(wm, n, :dhp, a)
    con_p = JuMP.@NLconstraint(wm.model, r * f_alpha(qp) - inv(L) * dhp <= 0.0)
    con(wm, n, :potential_loss_p)[a] = con_p
end

function objective_wf(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: StandardMICPForm
    JuMP.set_objective_sense(wm.model, MOI.FEASIBILITY_SENSE)
end

function objective_owf(wm::GenericWaterModel{T}) where T <: StandardMICPForm
    JuMP.set_objective_sense(wm.model, MOI.FEASIBILITY_SENSE)
end
