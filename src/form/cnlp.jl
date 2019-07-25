# Define CNLP (convex nonlinear programming) implementations of water distribution models.

function variable_head(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractCNLPForm
end

function variable_flow(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractCNLPForm
    variable_directed_flow(wm, n, bounded=false)
    variable_undirected_flow(wm, n, bounded=false)
end

function variable_pump(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractCNLPForm
end

function constraint_potential_loss_pipe(wm::GenericWaterModel{T}, n::Int, a::Int, alpha, f_id, t_id, len, r_min) where T <: AbstractCNLPForm
end

function constraint_head_difference(wm::GenericWaterModel{T}, n::Int, a::Int, f_id, t_id, head_fr, head_to) where T <: AbstractCNLPForm
end

function constraint_flow_direction_selection(wm::GenericWaterModel{T}, n::Int, a::Int) where T <: AbstractCNLPForm
end

function constraint_potential_loss_ub_pipe(wm::GenericWaterModel{T}, n::Int, a::Int, alpha, len, r_max) where T <: AbstractCNLPForm
end

function constraint_potential_loss_pump(wm::GenericWaterModel{T}, n::Int, a::Int, f_id::Int, t_id::Int) where T <: AbstractCNLPForm
end

function objective_wf(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: StandardCNLPForm
    linear_expr = JuMP.@expression(wm.model, 0.0)
    linear_expr_start = 0.0

    for (i, reservoir) in ref(wm, n, :reservoirs)
        for (a, link) in filter(a -> i == a.second["f_id"], ref(wm, n, :links))
            qp = var(wm, n, :qp, a)
            qn = var(wm, n, :qn, a)
            linear_expr -= reservoir["head"] * (qp - qn)
            linear_expr_start -= reservoir["head"] * (JuMP.start_value(qp) - JuMP.start_value(qn))
        end

        for (a, link) in filter(a -> i == a.second["t_id"], ref(wm, n, :links))
            qp = var(wm, n, :qp, a)
            qn = var(wm, n, :qn, a)
            linear_expr -= reservoir["head"] * (qn - qp)
            linear_expr_start -= reservoir["head"] * (JuMP.start_value(qn) - JuMP.start_value(qp))
        end
    end

    # TODO: Declare this variable somewhere else? Or even better, add something
    # to JuMP that allows for addition of affine and nonlinear expressions...
    linear_term = JuMP.@variable(wm.model, base_name="linear_objective_term", start=linear_expr_start)
    JuMP.@constraint(wm.model, linear_expr == linear_term)

    qn = var(wm, n, :qn)
    qp = var(wm, n, :qp)
    # Initialize the objective.
    objective_expr = JuMP.@NLexpression(wm.model, linear_term +
        sum(link["length"] * ref(wm, n, :resistance, a)[1] *
        (if_alpha(qn[a]) + if_alpha(qp[a]))
        for (a, link) in ref(wm, n, :links)))

    return JuMP.@NLobjective(wm.model, MOI.MIN_SENSE, objective_expr)
end
