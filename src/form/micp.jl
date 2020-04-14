# Define common MICP (mixed-integer convex program) implementations of water
# distribution constraints, which use directed flow variables.

function constraint_head_loss_pipe_des(wm::AbstractMICPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, pipe_resistances) 
    # Collect head difference variables.
    dhp, dhn = [var(wm, n, :dhp, a), var(wm, n, :dhn, a)]

    for (r_id, r) in enumerate(pipe_resistances)
        # Collect directed flow variables.
        qp, qn = [var(wm, n, :qp_des, a)[r_id], var(wm, n, :qn_des, a)[r_id]]

        # Build the relaxed head loss constraints.
        c_p = JuMP.@NLconstraint(wm.model, r * head_loss(qp) <= inv(L) * dhp)
        c_n = JuMP.@NLconstraint(wm.model, r * head_loss(qn) <= inv(L) * dhn)

        # Append the constraint array.
        append!(con(wm, n, :head_loss, a), [c_p, c_n])
    end
end

"Pump head gain constraint when the pump status is ambiguous."
function constraint_pump_head_gain(wm::AbstractMICPModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Gather flow and head gain variables.
    g = var(wm, n, :g, a)
    qp, x_pump = [var(wm, n, :qp, a), var(wm, n, :x_pump, a)]

    # Gather head-related variables and data.
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]
    dhp, dhn = [var(wm, n, :dhp, a), var(wm, n, :dhn, a)]
    dhp_ub, dhn_ub = [JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)]

    # Define the (relaxed) head gain caused by the pump.
    g_expr = curve_fun[1]*qp^2 + curve_fun[2]*qp + curve_fun[3]*x_pump
    c = JuMP.@constraint(wm.model, g_expr >= g) # Concavified.

    # Append the constraint array.
    append!(con(wm, n, :head_gain, a), [c])
end

function constraint_check_valve_head_loss(wm::AbstractMICPModel, n::Int, a::Int, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    # Gather flow- and check valve-related variables.
    qp, qn = [var(wm, n, :qp, a), var(wm, n, :qn, a)]
    x_cv = var(wm, n, :x_cv, a)

    # Gather head variables and upper bound data.
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]
    dhp, dhn = [var(wm, n, :dhp, a), var(wm, n, :dhn, a)]
    dhp_ub, dhn_ub = [JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)]

    # Add constraints for flow in the positive and negative directions.
    lhs = JuMP.@NLexpression(wm.model, r*head_loss(qp) - inv(L)*dhp)
    c_p = JuMP.@NLconstraint(wm.model, lhs <= dhp_ub * (1.0 - x_cv))
    c_n = JuMP.@NLconstraint(wm.model, dhn <= dhn_ub * (1.0 - x_cv))

    # Append the constraint array.
    append!(con(wm, n, :head_loss)[a], [c_p, c_n])
end

function constraint_head_loss_pipe(wm::AbstractMICPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    # Gather common variables.
    qp, qn = [var(wm, n, :qp, a), var(wm, n, :qn, a)]
    dhp, dhn = [var(wm, n, :dhp, a), var(wm, n, :dhn, a)]

    # Add constraints for head loss in the positive and negative directions.
    c_p = JuMP.@NLconstraint(wm.model, r * head_loss(qp) <= inv(L) * dhp)
    c_n = JuMP.@NLconstraint(wm.model, r * head_loss(qn) <= inv(L) * dhn)

    # Append the constraint array.
    append!(con(wm, n, :head_loss)[a], [c_p, c_n])
end

function objective_wf(wm::AbstractMICPModel)
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end

function objective_owf(wm::AbstractMICPModel)
    objective_wf(wm)
end
