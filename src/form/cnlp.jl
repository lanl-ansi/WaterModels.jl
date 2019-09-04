# Define CNLP (convex nonlinear programming) implementations of water distribution models.

function variable_head(wm::AbstractCNLPModel, n::Int=wm.cnw)
end

function variable_flow(wm::AbstractCNLPModel, n::Int=wm.cnw)
    variable_directed_flow(wm, n, bounded=false)
    variable_undirected_flow(wm, n)
end

function variable_pump(wm::AbstractCNLPModel, n::Int=wm.cnw)
end

function constraint_potential_loss_pipe(wm::AbstractCNLPModel, n::Int, a::Int, alpha, node_fr, node_to, len, r_min)
end

function constraint_head_difference(wm::AbstractCNLPModel, n::Int, a::Int, node_fr, node_to, head_fr, head_to)
end

function constraint_flow_direction_selection(wm::AbstractCNLPModel, n::Int, a::Int)
end

function constraint_potential_loss_ub_pipe(wm::AbstractCNLPModel, n::Int, a::Int, alpha, len, r_max)
end

function constraint_potential_loss_pump(wm::AbstractCNLPModel, n::Int, a::Int, node_fr::Int, node_to::Int)
end

function objective_wf(wm::AbstractCNLPModel, n::Int=wm.cnw)
    q_p = var(wm, n, :qp)
    q_n = var(wm, n, :qn)
    q_r = var(wm, n, :q_r)

    return JuMP.@NLobjective(wm.model, _MOI.MIN_SENSE,
        sum(-q_r[i] * res["head"] for (i, res) in ref(wm, n, :reservoirs)) +
        sum(link["length"] * ref(wm, n, :resistance, a)[1] *
        (head_loss(q_p[a]) + head_loss(q_n[a]))
        for (a, link) in ref(wm, n, :links)))
end
