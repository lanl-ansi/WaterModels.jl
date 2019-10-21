# Define MICP (mixed-integer convex program) implementations of water distribution constraints.

function constraint_head_loss_pipe_ne(wm::AbstractMICPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, pipe_resistances) 
    for (r_id, r) in enumerate(pipe_resistances)
        dhp = var(wm, n, :dhp, a)
        qp_ne = var(wm, n, :qp_ne, a)[r_id]
        cp = JuMP.@NLconstraint(wm.model, r * head_loss(qp_ne) <= inv(L) * dhp)

        dhn = var(wm, n, :dhn, a)
        qn_ne = var(wm, n, :qn_ne, a)[r_id]
        cn = JuMP.@NLconstraint(wm.model, r * head_loss(qn_ne) <= inv(L) * dhn)

        append!(con(wm, n, :head_loss, a), [cp, cn])
    end
end

"Pump head gain constraint when the pump status is ambiguous."
function constraint_head_gain_pump(wm::AbstractMICPModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Fix reverse flow variable to zero (since this is a pump).
    JuMP.fix(var(wm, n, :qn, a), 0.0, force=true)

    # Gather common variables.
    qp = var(wm, n, :qp, a)
    qp_ub = JuMP.has_upper_bound(qp) ? JuMP.upper_bound(qp) : 10.0

    g = var(wm, n, :g, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)
    x_pump = var(wm, n, :x_pump, a)

    # Define the (relaxed) head gain caused by the pump.
    con_1 = JuMP.@NLconstraint(wm.model, curve_fun[1]*qp^2 + curve_fun[2]*qp + curve_fun[3]*x_pump <= g)

    # If the pump is off, decouple the head difference relationship.
    con_2 = JuMP.@constraint(wm.model, h_j - h_i - g <= 1.0e6 * (1 - x_pump))
    con_3 = JuMP.@constraint(wm.model, h_j - h_i - g >= -1.0e6 * (1 - x_pump))

    # If the pump is off, the flow along the pump must be zero.
    con_4 = JuMP.@constraint(wm.model, qp <= qp_ub * x_pump)
    con_5 = JuMP.@constraint(wm.model, qp >= 1.0e-6 * x_pump)

    # Append the constraint array.
    con(wm, n, :head_gain)[a] = [con_1, con_2, con_3, con_4, con_5]
end

"Pump head gain constraint when the pump is forced to be on."
function constraint_head_gain_pump_on(wm::AbstractMICPModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Fix reverse flow variable to zero (since this is a pump).
    JuMP.fix(var(wm, n, :qn, a), 0.0, force=true)

    # Gather common variables.
    qp = var(wm, n, :qp, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)

    # Define the head difference relationship when the pump is on (h_j >= h_i).
    c = JuMP.@NLconstraint(wm.model, curve_fun[1]*qp^2 + curve_fun[2]*qp + curve_fun[3] <= (h_j - h_i))
    con(wm, n, :head_gain)[a] = [c]
end

function constraint_head_loss_check_valve(wm::AbstractMICPModel, n::Int, a::Int, node_fr::Int, node_to::Int, L::Float64, r::Float64) 
    # Gather common variables.
    qp = var(wm, n, :qp, a)
    qn = var(wm, n, :qn, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)
    x_cv = var(wm, n, :x_cv, a)

    lhs_p = JuMP.@NLexpression(wm.model, r * head_loss(qp) - inv(L) * (h_i - h_j))
    c_p = JuMP.@NLconstraint(wm.model, lhs_p <= 1.0e3 * (1.0 - x_cv))
    lhs_n = JuMP.@NLexpression(wm.model, r * head_loss(qn) - inv(L) * (h_j - h_i))
    c_n = JuMP.@NLconstraint(wm.model, lhs_n <= 1.0e3 * x_cv)

    append!(con(wm, n, :head_loss)[a], [c_p, c_n])
end

function constraint_head_loss_pipe(wm::AbstractMICPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    qp = var(wm, n, :qp, a)
    dhp = var(wm, n, :dhp, a)
    cp = JuMP.@NLconstraint(wm.model, r * head_loss(qp) <= inv(L) * dhp)

    qn = var(wm, n, :qn, a)
    dhn = var(wm, n, :dhn, a)
    cn = JuMP.@NLconstraint(wm.model, r * head_loss(qn) <= inv(L) * dhn)

    append!(con(wm, n, :head_loss)[a], [cp, cn])
end

function objective_wf(wm::AbstractMICPModel)
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end

function objective_owf(wm::AbstractMICPModel)
    objective_wf(wm)
end
