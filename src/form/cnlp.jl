# Define CNLP (convex nonlinear programming) implementations of water distribution models.

function variable_head(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractCNLPForm
end

function variable_flow(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractCNLPForm
    variable_directed_flow(wm, n, bounded=false)
    variable_undirected_flow(wm, n)
end

function variable_pump(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractCNLPForm
end

function constraint_potential_loss_pipe(wm::GenericWaterModel{T}, n::Int, a::Int, alpha, node_fr, node_to, len, r_min) where T <: AbstractCNLPForm
end

function constraint_head_difference(wm::GenericWaterModel{T}, n::Int, a::Int, node_fr, node_to, head_fr, head_to) where T <: AbstractCNLPForm
end

function constraint_flow_direction_selection(wm::GenericWaterModel{T}, n::Int, a::Int) where T <: AbstractCNLPForm
end

function constraint_potential_loss_ub_pipe(wm::GenericWaterModel{T}, n::Int, a::Int, alpha, len, r_max) where T <: AbstractCNLPForm
end

function constraint_potential_loss_pump(wm::GenericWaterModel{T}, n::Int, a::Int, node_fr::Int, node_to::Int) where T <: AbstractCNLPForm
end

function objective_wf(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: StandardCNLPForm
    q_p = var(wm, n, :qp)
    q_n = var(wm, n, :qn)
    q_r = var(wm, n, :q_r)

    return JuMP.@NLobjective(wm.model, MOI.MIN_SENSE,
        sum(-q_r[i] * res["head"] for (i, res) in ref(wm, n, :reservoirs)) +
        sum(link["length"] * ref(wm, n, :resistance, a)[1] *
        (head_loss(q_p[a]) + head_loss(q_n[a]))
        for (a, link) in ref(wm, n, :links)))
end
