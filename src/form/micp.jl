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

function constraint_potential_loss_pump(wm::GenericWaterModel{T}, n::Int, a::Int, f_id::Int, t_id::Int) where T <: AbstractMICPForm
end

function constraint_potential_loss_pipe_ne(wm::GenericWaterModel{T}, n::Int, a::Int, alpha, f_id, t_id, len, pipe_resistances) where T <: AbstractMICPForm
    if !haskey(con(wm, n), :potential_loss_n_ne)
        con(wm, n)[:potential_loss_n_ne] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        con(wm, n)[:potential_loss_p_ne] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    con(wm, n, :potential_loss_n_ne)[a] = Dict{Int, JuMP.ConstraintRef}()
    con(wm, n, :potential_loss_p_ne)[a] = Dict{Int, JuMP.ConstraintRef}()

    L = len

    #var(wm, n, :qn_ne, (a, r_id)) =
    for (r_id, r) in enumerate(pipe_resistances)
        qn_ne = var(wm, n, :qn_ne, a)[r_id]
        dhn = var(wm, n, :dhn, a)
        con_n = JuMP.@NLconstraint(wm.model, r * head_loss(qn_ne) - inv(L) * dhn <= 0.0)
        con(wm, n, :potential_loss_n_ne, a)[r_id] = con_n

        qp_ne = var(wm, n, :qp_ne, a)[r_id]
        dhp = var(wm, n, :dhp, a)
        con_p = JuMP.@NLconstraint(wm.model, r * head_loss(qp_ne) - inv(L) * dhp <= 0.0)
        con(wm, n, :potential_loss_p_ne, a)[r_id] = con_p
    end
end

function constraint_potential_loss_pipe(wm::GenericWaterModel{T}, n::Int, a::Int, alpha, f_id, t_id, len, r_min) where T <: AbstractMICPForm
    if !haskey(con(wm, n), :potential_loss_n)
        con(wm, n)[:potential_loss_n] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:potential_loss_p] = Dict{Int, JuMP.ConstraintRef}()
    end

    L = len
    r = r_min

    qn = var(wm, n, :qn, a)
    dhn = var(wm, n, :dhn, a)
    con_n = JuMP.@NLconstraint(wm.model, r * head_loss(qn) - inv(L) * dhn <= 0.0)
    con(wm, n, :potential_loss_n)[a] = con_n

    qp = var(wm, n, :qp, a)
    dhp = var(wm, n, :dhp, a)
    con_p = JuMP.@NLconstraint(wm.model, r * head_loss(qp) - inv(L) * dhp <= 0.0)
    con(wm, n, :potential_loss_p)[a] = con_p
end

function objective_wf(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: StandardMICPForm
    JuMP.set_objective_sense(wm.model, MOI.FEASIBILITY_SENSE)
end

function objective_owf(wm::GenericWaterModel{T}) where T <: StandardMICPForm
    JuMP.set_objective_sense(wm.model, MOI.FEASIBILITY_SENSE)
end
