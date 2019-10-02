# Define CNLP (convex nonlinear programming) implementations of water network models.

"Head variables do not exist within the CNLP model."
function variable_head(wm::AbstractCNLPModel, n::Int=wm.cnw; bounded::Bool=true)
    if length(ids(wm, n, :tank)) >= 1
        Memento.error(_LOGGER, "AbstractCNLPModel does not support tank components.")
    end
end

"Pump variables do not exist within the CNLP model."
function variable_pump(wm::AbstractCNLPModel, n::Int=wm.cnw) end

"The source head constraint does not exist within the CNLP model."
function constraint_source_head(wm::AbstractCNLPModel, n::Int, i::Int, h_src::Float64) end

"The head loss constraint does not exist within the CNLP model."
function constraint_head_loss_pipe(wm::AbstractCNLPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, len::Float64, r::Float64) end

"The head difference constraint does not exist within the CNLP model."
function constraint_head_difference(wm::AbstractCNLPModel, n::Int, a::Int, node_fr::Int, node_to::Int, head_fr, head_to) end

"The flow direction selection constraint does not exist within the CNLP model."
function constraint_flow_direction_selection(wm::AbstractCNLPModel, n::Int, a::Int) end

"The head loss upper bound constraint does not exist within the CNLP model."
function constraint_head_loss_ub_pipe(wm::AbstractCNLPModel, n::Int, a::Int, alpha::Float64, len::Float64, r::Float64) end

"The pump head loss constraint does not exist within the CNLP model."
function constraint_head_loss_pump(wm::AbstractCNLPModel, n::Int, a::Int, node_fr::Int, node_to::Int) end

"The objective function of the CNLP model results in first-order optimality
conditions that are similar to the full, nonconvex head loss constraints."
function objective_wf(wm::AbstractCNLPModel, n::Int=wm.cnw)
    # TODO: How can the below result in pump and tank KKT conditions?
    return JuMP.@NLobjective(wm.model, _MOI.MIN_SENSE,
        sum(sum(-wm.var[:nw][n][:qr][i] * res["head"] for (i, res) in
        ref(wm, n, :reservoir)) + sum(pipe["length"] *
        ref(wm, n, :resistance, a)[1] * (head_loss(wm.var[:nw][n][:qp][a]) +
        head_loss(wm.var[:nw][n][:qn][a])) for (a, pipe) in ref(wm, n, :pipe))
        for n in nw_ids(wm)))
end

function objective_cwf(wm::AbstractCNLPModel, n::Int=wm.cnw)
    return objective_wf(wm, n)
end

"Pump head gain constraint when the pump status is ambiguous."
function constraint_head_gain_pump(wm::AbstractCNLPModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Fix reverse flow variable to zero (since this is a pump).
    qn = var(wm, n, :qn, a)
    JuMP.fix(qn, 0.0, force=true)

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
function constraint_head_gain_pump_on(wm::AbstractCNLPModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Fix reverse flow variable to zero (since this is a pump).
    qn = var(wm, n, :qn, a)
    JuMP.fix(qn, 0.0, force=true)

    # Gather common variables.
    qp = var(wm, n, :qp, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)

    # Define the head difference relationship when the pump is on (h_j >= h_i).
    con_1 = JuMP.@NLconstraint(wm.model, curve_fun[1]*qp^2 + curve_fun[2]*qp + curve_fun[3] <= (h_j - h_i))
    con(wm, n, :head_gain)[a] = [con_1]
end
