# Define MICP (mixed-integer convex program) implementations of water distribution constraints.

function constraint_head_loss_pipe_ne(wm::AbstractMICPModel, n::Int, a::Int,
    alpha::Float64, node_fr::Int, node_to::Int, L::Float64, pipe_resistances) 
    for (r_id, r) in enumerate(pipe_resistances)
        dhp = var(wm, n, :dhp, a)
        qp = var(wm, n, :qp_ne, a)[r_id]
        cp = JuMP.@NLconstraint(wm.model, r * head_loss(qp) <= inv(L) * dhp)

        dhn = var(wm, n, :dhn, a)
        qn = var(wm, n, :qn_ne, a)[r_id]
        cn = JuMP.@NLconstraint(wm.model, r * head_loss(qn) <= inv(L) * dhn)

        append!(con(wm, n, :head_loss, a), [cp, cn])
    end
end

"Pump head gain constraint when the pump status is ambiguous."
function constraint_head_gain_pump(wm::AbstractMICPModel, n::Int, a::Int,
    node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Gather flow-related variables and data.
    qp = var(wm, n, :qp, a)
    qp_ub = JuMP.upper_bound(qp)

    # Gather pump-related variables.
    g = var(wm, n, :g, a)
    x_p = var(wm, n, :x_pump, a)

    # Gather head variables.
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)

    # Gather head difference-related variables and data.
    dhp = var(wm, n, :dhp, a)
    dhp_ub = JuMP.upper_bound(dhp)
    dhn = var(wm, n, :dhn, a)
    dhn_ub = JuMP.upper_bound(dhn)

    # Define the (relaxed) head gain caused by the pump.
    g_expr = curve_fun[1]*qp^2 + curve_fun[2]*qp + curve_fun[3]*x_p
    c_1 = JuMP.@constraint(wm.model, g_expr <= g) # Convexified.

    # If the pump is off, decouple the head difference relationship.
    c_2 = JuMP.@constraint(wm.model, dhn <= g + dhn_ub * (1.0 - x_p))
    c_3 = JuMP.@constraint(wm.model, dhp <= dhp_ub * (1.0 - x_p))

    # If the pump is off, the flow along the pump must be zero.
    c_4 = JuMP.@constraint(wm.model, qp <= qp_ub * x_p)
    c_5 = JuMP.@constraint(wm.model, qp >= 6.31465679e-6 * x_p)

    # Append the constraint array.
    con(wm, n, :head_gain)[a] = [c_1, c_2, c_3, c_4, c_5]
end

"Pump head gain constraint when the pump is forced to be on."
function constraint_head_gain_pump_on(wm::AbstractMICPModel, n::Int, a::Int,
    node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Gather common variables.
    qp = var(wm, n, :qp, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)

    # Define the head difference relationship when the pump is on (h_j >= h_i).
    g_expr = curve_fun[1]*qp^2 + curve_fun[2]*qp + curve_fun[3]*x_p
    c = JuMP.@NLconstraint(wm.model, g_expr <= h_j - h_i)

    con(wm, n, :head_gain)[a] = [c]
end

function constraint_head_loss_check_valve(wm::AbstractMICPModel, n::Int, a::Int,
    node_fr::Int, node_to::Int, L::Float64, r::Float64) 
    # Gather flow-related variables.
    qp = var(wm, n, :qp, a)
    qn = var(wm, n, :qn, a)
    x_cv = var(wm, n, :x_cv, a)

    # Gather head variables.
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)

    # Gather head difference-related variables and data.
    dhp = var(wm, n, :dhp, a)
    dhp_ub = JuMP.upper_bound(dhp)
    dhn = var(wm, n, :dhn, a)
    dhn_ub = JuMP.upper_bound(dhn)

    # Add constraints for flow in the positive and negative directions.
    lhs_p = JuMP.@NLexpression(wm.model, r*head_loss(qp) - inv(L)*dhp)
    c_p = JuMP.@NLconstraint(wm.model, lhs_p <= dhp_ub * (1.0 - x_cv))
    c_n = JuMP.@NLconstraint(wm.model, dhn <= dhn_ub * (1.0 - x_cv))
    #lhs_n = JuMP.@NLexpression(wm.model, r*head_loss(qn) - inv(L)*dhn)

    append!(con(wm, n, :head_loss)[a], [c_p, c_n])
end

function constraint_head_loss_pipe(wm::AbstractMICPModel, n::Int, a::Int,
    alpha::Float64, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    qp = var(wm, n, :qp, a)
    dhp = var(wm, n, :dhp, a)
    cp = JuMP.@NLconstraint(wm.model, r*head_loss(qp) <= inv(L)*dhp)

    qn = var(wm, n, :qn, a)
    dhn = var(wm, n, :dhn, a)
    cn = JuMP.@NLconstraint(wm.model, r*head_loss(qn) <= inv(L)*dhn)

    append!(con(wm, n, :head_loss)[a], [cp, cn])
end

function objective_wf(wm::AbstractMICPModel)
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end

function objective_owf(wm::AbstractMICPModel)
    objective_wf(wm)
end
