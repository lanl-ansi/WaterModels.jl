# Constraints and variables common to all formulations with directed flows. In
# these formulations, the variables qp correspond to flow from i to j, and the
# variables qn correspond to flow from j to i. That is, when qp is nonzero, qn
# should be zero, and when qn is nonzero, qp should be zero.

"Create head-related variables for formulations that use binary direction
variables. These head variables include total hydraulic head variables as well
as positive (i to j) and negative (j to i) head difference variables."
function variable_head(wm::AbstractIntegerDirectedFlowForms; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    # Initialize variables associated with heads.
    h = var(wm, nw)[:h] = JuMP.@variable(wm.model,
        [i in ids(wm, nw, :node)], base_name="$(nw)_h",
        start=comp_start_value(ref(wm, nw, :node, i), "h_start"))

    # Initialize variables associated with positive head differences.
    dhp = var(wm, nw)[:dhp] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :link)], lower_bound=0.0, base_name="$(nw)_dhp",
        start=comp_start_value(ref(wm, nw, :link, a), "dhp_start"))

    # Initialize variables associated with negative head differences.
    dhn = var(wm, nw)[:dhn] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :link)], lower_bound=0.0, base_name="$(nw)_dhn",
        start=comp_start_value(ref(wm, nw, :link, a), "dhn_start"))

    if bounded # Bound head-related variables if required.
        # Get the head bound variables.
        h_lb, h_ub = calc_head_bounds(wm, nw)

        # Set lower and upper bounds on heads.
        for (i, node) in ref(wm, nw, :node)
            JuMP.set_lower_bound(h[i], h_lb[i])
            JuMP.set_upper_bound(h[i], h_ub[i])
        end

        # Set lower and upper bounds on head differences.
        for (a, link) in ref(wm, nw, :link)
            i, j = [link["node_fr"], link["node_to"]]
            JuMP.set_upper_bound(dhp[a], max(0.0, h_ub[i] - h_lb[j]))
            JuMP.set_upper_bound(dhn[a], max(0.0, h_ub[j] - h_lb[i]))
        end
    end

    # Report back head values as part of the solution.
    report && sol_component_value(wm, nw, :node, :h, ids(wm, nw, :node), h)
end

"Create flow variables for CNLP model."
function variable_flow(wm::AbstractContinuousDirectedFlowForms; nw::Int=wm.cnw, bounded::Bool=true)
    variable_flow_common(wm, nw=nw, bounded=bounded)
end

"Create flow variables for formulations with binary direction variables."
function variable_flow(wm::AbstractIntegerDirectedFlowForms; nw::Int=wm.cnw, bounded::Bool=true)
    variable_flow_common(wm, nw=nw, bounded=bounded)
    variable_flow_direction(wm, nw=nw)
end

"Create flow variables that are common to all directed flow models."
function variable_flow_common(wm::AbstractDirectedFlowModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    # Initialize variables associated with positive flows.
    qp = var(wm, nw)[:qp] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :link_fixed)], lower_bound=0.0, base_name="$(nw)_qp",
        start=comp_start_value(ref(wm, nw, :link_fixed, a), "qp_start"))

    # Initialize variables associated with negative flows.
    qn = var(wm, nw)[:qn] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :link_fixed)], lower_bound=0.0, base_name="$(nw)_qn",
        start=comp_start_value(ref(wm, nw, :link_fixed, a), "qn_start"))

    if bounded # Bound flow-related variables if desired.
        q_lb, q_ub = calc_flow_bounds(wm, nw)

        for (a, link) in ref(wm, nw, :link_fixed)
            JuMP.set_upper_bound(qp[a], max(0.0, maximum(q_ub[a])))
            JuMP.set_upper_bound(qn[a], max(0.0, -minimum(q_lb[a])))
        end
    end

    # Create expressions capturing the relationships among q, qp, and qn.
    q = var(wm, nw)[:q] = JuMP.@expression(wm.model,
        [a in ids(wm, nw, :link_fixed)],
        var(wm, nw, :qp, a) - var(wm, nw, :qn, a))

    # Report back flow values as part of the solution.
    report && sol_component_value(wm, nw, :link, :q, ids(wm, nw, :link_fixed), q)
end

"Initialize variables associated with flow direction. If this variable is equal
to one, the flow direction is from i to j. If it is equal to zero, the flow
direction is from j to i."
function variable_flow_direction(wm::AbstractDirectedFlowModel; nw::Int=wm.cnw)
    var(wm, nw)[:x_dir] = JuMP.@variable(wm.model, [a in ids(wm, nw, :link)],
        base_name="$(nw)_x_dir", binary=true,
        start=comp_start_value(ref(wm, nw, :link, a), "x_dir_start"))
end

"Create network design flow variables for directed flow formulations."
function variable_flow_des(wm::AbstractDirectedFlowModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    # Create dictionary for undirected design flow variables (qp_des and qn_des).
    qp_des = var(wm, nw)[:qp_des] = Dict{Int,Array{JuMP.VariableRef}}()
    qn_des = var(wm, nw)[:qn_des] = Dict{Int,Array{JuMP.VariableRef}}()

    # Initialize the variables. (The default start value of 1.0e-6 is crucial.)
    for a in ids(wm, nw, :link_des)
        var(wm, nw, :qp_des)[a] = JuMP.@variable(wm.model,
            [r in 1:length(ref(wm, nw, :resistance, a))], lower_bound=0.0,
            base_name="$(nw)_qp_des[$(a)]",
            start=comp_start_value(ref(wm, nw, :link_des, a), "qp_des_start", r, 1.0e-6))

        var(wm, nw, :qn_des)[a] = JuMP.@variable(wm.model,
            [r in 1:length(ref(wm, nw, :resistance, a))], lower_bound=0.0,
            base_name="$(nw)_qn_des[$(a)]",
            start=comp_start_value(ref(wm, nw, :link_des, a), "qn_des_start", r, 1.0e-6))
    end

    if bounded # If the variables are bounded, apply the bounds.
        q_lb, q_ub = calc_flow_bounds(wm, nw)

        for a in ids(wm, nw, :link_des)
            for r in 1:length(ref(wm, nw, :resistance, a))
                JuMP.set_upper_bound(qp_des[a][r], max(0.0, q_ub[a][r]))
                JuMP.set_upper_bound(qn_des[a][r], max(0.0, -q_lb[a][r]))
            end
        end
    end

    # Create expressions capturing the relationships among q, qp_des, and qn_des.
    q = var(wm, nw)[:q] = JuMP.@expression(
        wm.model, [a in ids(wm, nw, :link_des)],
        sum(var(wm, nw, :qp_des, a)) - sum(var(wm, nw, :qn_des, a)))

    # Initialize the solution reporting data structures.
    report && sol_component_value(wm, nw, :link, :q, ids(wm, nw, :link_des), q)

    # Create resistance binary variables.
    variable_resistance(wm, nw=nw)
end

function constraint_check_valve_common(wm::AbstractDirectedFlowModel, n::Int, a::Int, node_fr::Int, node_to::Int, head_fr, head_to)
    # Get common flow variables.
    qp, x_cv = [var(wm, n, :qp, a), var(wm, n, :x_cv, a)]

    # If the check valve is open, flow must be appreciably nonnegative.
    c_1 = JuMP.@constraint(wm.model, qp <= JuMP.upper_bound(qp) * x_cv)
    c_2 = JuMP.@constraint(wm.model, qp >= 6.31465679e-6 * x_cv)

    # Get common head variables and associated data.
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]
    dhp, dhn = [var(wm, n, :dhp, a), var(wm, n, :dhn, a)]
    dhp_ub, dhn_ub = [JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)]

    # When the check valve is open, negative head loss is not possible.
    c_3 = JuMP.@constraint(wm.model, dhn <= dhn_ub * (1.0 - x_cv))

    # When the check valve is closed, positive head loss is not possible.
    c_4 = JuMP.@constraint(wm.model, dhp <= dhp_ub * x_cv)

    # Append the constraint array.
    append!(con(wm, n, :check_valve, a), [c_1, c_2, c_3, c_4])
end

function constraint_pipe_common(wm::AbstractDirectedFlowModel, n::Int, a::Int, node_fr::Int, node_to::Int, head_fr, head_to)
    # Get common flow variables and associated data.
    x_dir = var(wm, n, :x_dir, a)
    qp, qn = [var(wm, n, :qp, a), var(wm, n, :qn, a)]
    qp_ub, qn_ub = [JuMP.upper_bound(qp), JuMP.upper_bound(qn)]

    # Constrain directed flow variables based on direction.
    c_1 = JuMP.@constraint(wm.model, qp <= qp_ub * x_dir)
    c_2 = JuMP.@constraint(wm.model, qn <= qn_ub * (1.0 - x_dir))

    # Get common head variables and associated data.
    dhp, dhn = [var(wm, n, :dhp, a), var(wm, n, :dhn, a)]
    h_i = head_fr == nothing ? var(wm, n, :h, node_fr) : head_fr
    h_j = head_to == nothing ? var(wm, n, :h, node_to) : head_to
    dhp_ub, dhn_ub = [JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)]

    # Constraint directed head variables based on direction.
    c_3 = JuMP.@constraint(wm.model, dhp <= dhp_ub * x_dir)
    c_4 = JuMP.@constraint(wm.model, dhn <= dhn_ub * (1.0 - x_dir))
    c_5 = JuMP.@constraint(wm.model, dhp - dhn == h_i - h_j)

    # Append the constraint array.
    append!(con(wm, n, :pipe)[a], [c_1, c_2, c_3, c_4, c_5])
end

function constraint_pump_common(wm::AbstractDirectedFlowModel, n::Int, a::Int, node_fr::Int, node_to::Int, head_fr, head_to, pc::Array{Float64})
    # Gather variables common to all formulations involving pumps.
    x_dir, x_pump = [var(wm, n, :x_dir, a), var(wm, n, :x_pump, a)]
    qp, g = [var(wm, n, :qp, a), var(wm, n, :g, a)]
    dhp, dhn = [var(wm, n, :dhp, a), var(wm, n, :dhn, a)]

    # If the pump is off, flow across the pump must be zero.
    qp_ub = JuMP.upper_bound(qp)
    c_1 = JuMP.@constraint(wm.model, qp <= qp_ub * x_pump)
    c_2 = JuMP.@constraint(wm.model, qp >= 6.31465679e-6 * x_pump)

    # If the pump is off, decouple the head difference relationship.
    dhp_ub, dhn_ub = [JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)]
    c_3 = JuMP.@constraint(wm.model, dhp <= dhp_ub * (1.0 - x_pump))
    c_4 = JuMP.@constraint(wm.model, dhn <= g + dhn_ub * (1.0 - x_pump))

    # If the pump is on, flow across the pump must be nonnegative.
    c_5 = JuMP.@constraint(wm.model, x_dir == x_pump)

    # Equate head difference variables with heads.
    h_i = head_fr == nothing ? var(wm, n, :h, node_fr) : head_fr
    h_j = head_to == nothing ? var(wm, n, :h, node_to) : head_to
    c_6 = JuMP.@constraint(wm.model, dhp - dhn == h_i - h_j)

    # Add a linear lower bound on the head gain approximation.
    g_1, g_2 = [pc[3], pc[1]*qp_ub^2 + pc[2]*qp_ub + pc[3]]
    g_lb_line = (g_2 - g_1) * inv(qp_ub) * qp + g_1 * x_pump
    c_7 = JuMP.@constraint(wm.model, g_lb_line <= g)

    # Append the constraint array.
    append!(con(wm, n, :pump, a), [c_1, c_2, c_3, c_4, c_5, c_6, c_7])
end

function constraint_flow_direction_selection_des(wm::AbstractDirectedFlowModel, n::Int, a::Int, pipe_resistances)
    x_dir = var(wm, n, :x_dir, a)

    for r_id in 1:length(pipe_resistances)
        qp, qn = [var(wm, n, :qp_des, a)[r_id], var(wm, n, :qn_des, a)[r_id]]
        c_p = JuMP.@constraint(wm.model, qp <= JuMP.upper_bound(qp) * x_dir)
        c_n = JuMP.@constraint(wm.model, qn <= JuMP.upper_bound(qn) * (1.0 - x_dir))
        append!(con(wm, n, :head_loss)[a], [c_p, c_n])
    end
end

function constraint_head_loss_ub_cv(wm::AbstractDirectedFlowModel, n::Int, a::Int, alpha::Float64, L::Float64, r::Float64)
    qp, dhp = [var(wm, n, :qp, a), var(wm, n, :dhp, a)]
    rhs_p = r * JuMP.upper_bound(qp)^(alpha - 1.0) * qp
    c_p = JuMP.@constraint(wm.model, inv(L) * dhp <= rhs_p)
    append!(con(wm, n, :head_loss)[a], [c_p])
end

function constraint_head_loss_ub_pipe(wm::AbstractDirectedFlowModel, n::Int, a::Int, alpha::Float64, L::Float64, r::Float64)
    qp, dhp = [var(wm, n, :qp, a), var(wm, n, :dhp, a)]
    rhs_p = r * JuMP.upper_bound(qp)^(alpha - 1.0) * qp
    c_p = JuMP.@constraint(wm.model, inv(L) * dhp <= rhs_p)

    qn, dhn = [var(wm, n, :qn, a), var(wm, n, :dhn, a)]
    rhs_n = r * JuMP.upper_bound(qn)^(alpha - 1.0) * qn
    c_n = JuMP.@constraint(wm.model, inv(L) * dhn <= rhs_n)

    append!(con(wm, n, :head_loss)[a], [c_p, c_n])
end

function constraint_head_loss_ub_pipe_des(wm::AbstractDirectedFlowModel, n::Int, a::Int, alpha::Float64, L::Float64, pipe_resistances)
    dhp = var(wm, n, :dhp, a)
    qp_des = var(wm, n, :qp_des, a)
    qp_des_ub = JuMP.upper_bound.(qp_des)
    slopes_p = pipe_resistances .* qp_des_ub.^(alpha - 1.0)
    c_p = JuMP.@constraint(wm.model, inv(L)*dhp <= sum(slopes_p .* qp_des))

    dhn = var(wm, n, :dhn, a)
    qn_des = var(wm, n, :qn_des, a)
    qn_des_ub = JuMP.upper_bound.(qn_des)
    slopes_n = pipe_resistances .* qn_des_ub.^(alpha - 1.0)
    c_n = JuMP.@constraint(wm.model, inv(L)*dhn <= sum(slopes_n .* qn_des))

    append!(con(wm, n, :head_loss)[a], [c_p, c_n])
end

function constraint_resistance_selection_des(wm::AbstractDirectedFlowModel, n::Int, a::Int, pipe_resistances)
    c = JuMP.@constraint(wm.model, sum(var(wm, n, :x_res, a)) == 1.0)
    append!(con(wm, n, :head_loss)[a], [c])

    for (r_id, r) in enumerate(pipe_resistances)
        x_res = var(wm, n, :x_res, a)[r_id]

        qp_des = var(wm, n, :qp_des, a)[r_id]
        qp_ub = JuMP.upper_bound(qp_des)
        c_p = JuMP.@constraint(wm.model, qp_des <= qp_ub * x_res)

        qn_des = var(wm, n, :qn_des, a)[r_id]
        qn_ub = JuMP.upper_bound(qn_des)
        c_n = JuMP.@constraint(wm.model, qn_des <= qn_ub * x_res)

        append!(con(wm, n, :head_loss)[a], [c_p, c_n])
    end
end

"Constraint to ensure at least one direction is set to take flow away from a source."
function constraint_source_flow(wm::AbstractIntegerDirectedFlowForms, n::Int, i::Int, a_fr::Array{Tuple{Int,Int,Int}}, a_to::Array{Tuple{Int,Int,Int}})
    # Collect the required variables.
    x_dir = var(wm, n, :x_dir)

    # Get sums of inward- and outward-pointing direction variables.
    x_out = Array{JuMP.VariableRef}([x_dir[a] for (a, f, t) in a_fr])
    x_in = Array{JuMP.VariableRef}([x_dir[a] for (a, f, t) in a_to])

    # Add the source flow direction constraint.
    c = JuMP.@constraint(wm.model, sum(x_out) - sum(x_in) >= 1.0 - length(x_in))
    con(wm, n, :source_flow)[i] = c
end

"Constraint to ensure at least one direction is set to take flow to a junction with demand."
function constraint_sink_flow(wm::AbstractIntegerDirectedFlowForms, n::Int, i::Int, a_fr::Array{Tuple{Int,Int,Int}}, a_to::Array{Tuple{Int,Int,Int}})
    # Collect the required variables.
    x_dir = var(wm, n, :x_dir)

    # Get sums of inward- and outward-pointing direction variables.
    x_out = Array{JuMP.VariableRef}([x_dir[a] for (a, f, t) in a_fr])
    x_in = Array{JuMP.VariableRef}([x_dir[a] for (a, f, t) in a_to])

    # Add the sink flow direction constraint.
    c = JuMP.@constraint(wm.model, sum(x_in) - sum(x_out) >= 1.0 - length(x_out))
    con(wm, n, :sink_flow)[i] = c
end

function constraint_energy_conservation(wm::AbstractDirectedFlowModel, n::Int, r, L, alpha)
end
