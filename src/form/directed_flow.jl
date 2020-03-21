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
        alpha = ref(wm, nw, :alpha)
        q_lb, q_ub = calc_flow_bounds(wm, nw)

        for (a, link) in ref(wm, nw, :link_fixed)
            JuMP.set_upper_bound(qp[a], max(0.0, maximum(q_ub[a])))
            JuMP.set_upper_bound(qn[a], max(0.0, -minimum(q_lb[a])))
        end
    end

    q = var(wm, nw)[:q] = JuMP.@expression(wm.model,
        [a in ids(wm, nw, :link_fixed)],
        var(wm, nw, :qp, a) - var(wm, nw, :qn, a))

    report && sol_component_value(wm, nw, :link, :q,
        ids(wm, nw, :link_fixed), q)
end

"Initialize variables associated with flow direction. If this variable is equal
to one, the flow direction is from i to j. If it is equal to zero, the flow
direction is from j to i."
function variable_flow_direction(wm::AbstractDirectedFlowModel; nw::Int=wm.cnw)
    var(wm, nw)[:x_dir] = JuMP.@variable(wm.model, [a in ids(wm, nw, :link)],
        base_name="$(nw)_x_dir", binary=true,
        start=comp_start_value(ref(wm, nw, :link, a), "x_dir_start"))
end

"Create network expansion flow variables for directed flow formulations."
function variable_flow_ne(wm::AbstractIntegerDirectedFlowForms; nw::Int=wm.cnw, bounded::Bool=true)
    # Create directed qp_ne and qn_ne variables.
    bounded ? variable_flow_bounded_ne(wm, nw=nw) :
        variable_flow_unbounded_ne(wm, nw=nw)

    # Create expressions capturing the relationships among q, qp, and qn.
    var(wm, nw)[:q] = JuMP.@expression(wm.model, [a in ids(wm, nw, :link_ne)],
        sum(var(wm, nw, :qp_ne, a)) - sum(var(wm, nw, :qn_ne, a)))

    # Create resistance binary variables.
    variable_resistance(wm, nw=nw)
end

"Create network expansion flow variables for the CNLP formulation."
function variable_flow_ne(wm::AbstractContinuousDirectedFlowForms; nw::Int=wm.cnw, bounded::Bool=true)
    if length(ref(wm, nw, :link_ne)) > 0
        Memento.error(_LOGGER, "AbstractCNLPModel does not support network expansion formulations.")
    end
end

"Create network expansion bounded flow variables for directed flow formulations."
function variable_flow_bounded_ne(wm::AbstractDirectedFlowModel; nw::Int=wm.cnw)
    alpha = ref(wm, nw, :alpha)
    qn_ub, qp_ub = calc_directed_flow_upper_bounds(wm, alpha, nw)
    link_ids = ids(wm, nw, :link_ne)

    var(wm, nw)[:qp_ne] = Dict{Int,Array{JuMP.VariableRef}}(a => [] for a in link_ids)
    var(wm, nw)[:qn_ne] = Dict{Int,Array{JuMP.VariableRef}}(a => [] for a in link_ids)

    for a in link_ids
        n_r = length(ref(wm, nw, :resistance, a)) # Number of resistances.

        var(wm, nw, :qp_ne)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            lower_bound=0.0, upper_bound=qp_ub[a][r], base_name="qp_ne[$(nw)][$(a)]",
            start=get_start(ref(wm, nw, :link_ne), a, r, "qp_ne_start", 1.0e-6))

        var(wm, nw, :qn_ne)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            lower_bound=0.0, upper_bound=qn_ub[a][r], base_name="qn_ne[$(nw)][$(a)]",
            start=get_start(ref(wm, nw, :link_ne), a, r, "qn_ne_start", 1.0e-6))
    end
end

"Create network expansion unbounded flow variables for directed flow formulations."
function variable_flow_unbounded_ne(wm::AbstractDirectedFlowModel; nw::Int=wm.cnw)
    var(wm, nw)[:qp_ne] = Dict{Int,Array{JuMP.VariableRef}}(a => [] for a in ids(wm, nw, :link_ne))
    var(wm, nw)[:qn_ne] = Dict{Int,Array{JuMP.VariableRef}}(a => [] for a in ids(wm, nw, :link_ne))

    for a in ids(wm, nw, :link_ne)
        n_r = length(ref(wm, nw, :resistance, a)) # Number of resistances.

        var(wm, nw, :qp_ne)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            lower_bound=0.0, base_name="qp_ne[$(nw)][$(a)]",
            start=get_start(ref(wm, nw, :link_ne), a, r, "qp_ne_start", 1.0e-6))

        var(wm, nw, :qn_ne)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            lower_bound=0.0, base_name="qn_ne[$(nw)][$(a)]",
            start=get_start(ref(wm, nw, :link_ne), a, r, "qn_ne_start", 1.0e-6))
    end
end

function constraint_check_valve(wm::AbstractDirectedFlowModel, n::Int, a::Int, node_fr::Int, node_to::Int)
    # Fix negative flow variable to zero.
    JuMP.fix(var(wm, n, :qn, a), 0.0, force=true)

    # Gather common variables.
    x_cv = var(wm, n, :x_cv, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)

    # If the check valve is open, flow must be appreciably nonnegative.
    qp = var(wm, n, :qp, a)
    c_1 = JuMP.@constraint(wm.model, qp <= JuMP.upper_bound(qp) * x_cv)
    c_2 = JuMP.@constraint(wm.model, qp >= 6.31465679e-6 * x_cv)

    dhp = var(wm, n, :dhp, a)
    dhp_ub = JuMP.upper_bound(dhp)
    c_3 = JuMP.@constraint(wm.model, dhp <= x_cv * dhp_ub)

    dhn = var(wm, n, :dhn, a)
    dhn_ub = JuMP.upper_bound(dhn)
    c_4 = JuMP.@constraint(wm.model, dhn <= (1.0 - x_cv) * dhn_ub)

    c_e = JuMP.@constraint(wm.model, dhp - dhn == h_i - h_j)

    # TODO: These constraints seem to result in infeasibility in multiperiod Richmond case.
    #dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    #dh_ub = JuMP.upper_bound(h_i) - JuMP.lower_bound(h_j)
    #c_3 = JuMP.@constraint(wm.model, h_i - h_j >= (1.0 - x_cv) * dh_lb)
    #c_4 = JuMP.@constraint(wm.model, h_i - h_j <= x_cv * dh_ub)

    append!(con(wm, n, :check_valve)[a], [c_1, c_2, c_3, c_4, c_e])
end

"Create constraints for directed flow variables, based on flow direction
variables."
function constraint_flow_direction_selection(wm::AbstractDirectedFlowModel, n::Int, a::Int)
    x_dir = var(wm, n, :x_dir, a)

    qp = var(wm, n, :qp, a)
    qp_ub = JuMP.upper_bound(qp)
    cp = JuMP.@constraint(wm.model, qp <= qp_ub * x_dir)

    qn = var(wm, n, :qn, a)
    qn_ub = JuMP.upper_bound(qn)
    cn = JuMP.@constraint(wm.model, qn <= qn_ub * (1.0 - x_dir))

    append!(con(wm, n, :head_loss)[a], [cp, cn])
end

function constraint_flow_direction_selection_ne(wm::AbstractDirectedFlowModel, n::Int, a::Int, pipe_resistances)
    x_dir = var(wm, n, :x_dir, a)

    for r_id in 1:length(pipe_resistances)
        qp_ne = var(wm, n, :qp_ne, a)[r_id]
        qp_ne_ub = JuMP.has_upper_bound(qp_ne) ? JuMP.upper_bound(qp_ne) : 10.0
        cp = JuMP.@constraint(wm.model, qp_ne <= qp_ne_ub * x_dir)

        qn_ne = var(wm, n, :qn_ne, a)[r_id]
        qn_ne_ub = JuMP.has_upper_bound(qn_ne) ? JuMP.upper_bound(qn_ne) : 10.0
        cn = JuMP.@constraint(wm.model, qn_ne <= qn_ne_ub * (1.0 - x_dir))

        append!(con(wm, n, :head_loss)[a], [cp, cn])
    end
end

function constraint_head_difference(wm::AbstractDirectedFlowModel,
    n::Int, a::Int, node_fr::Int, node_to::Int, head_fr, head_to)
    x_dir = var(wm, n, :x_dir, a)

    dhp = var(wm, n, :dhp, a)
    dhp_ub = JuMP.upper_bound(dhp)
    cp = JuMP.@constraint(wm.model, dhp <= dhp_ub * x_dir)

    dhn = var(wm, n, :dhn, a)
    dhn_ub = JuMP.upper_bound(dhn)
    cn = JuMP.@constraint(wm.model, dhn <= dhn_ub * (1.0 - x_dir))

    h_i = head_fr == nothing ? var(wm, n, :h, node_fr) : head_fr
    h_j = head_to == nothing ? var(wm, n, :h, node_to) : head_to

    ce = JuMP.@constraint(wm.model, dhp - dhn == h_i - h_j)
    append!(con(wm, n, :head_loss)[a], [cp, cn, ce])
end

function constraint_head_loss_ub_pipe(wm::AbstractDirectedFlowModel,
    n::Int, a::Int, alpha::Float64, L::Float64, r::Float64)
    qp = var(wm, n, :qp, a)
    dhp = var(wm, n, :dhp, a)
    rhs_p = r*JuMP.upper_bound(qp)^(alpha - 1.0) * qp
    c_p = JuMP.@constraint(wm.model, inv(L)*dhp <= rhs_p)

    qn = var(wm, n, :qn, a)
    dhn = var(wm, n, :dhn, a)
    rhs_n = r*JuMP.upper_bound(qn)^(alpha - 1.0) * qn
    c_n = JuMP.@constraint(wm.model, inv(L)*dhn <= rhs_n)

    append!(con(wm, n, :head_loss)[a], [c_p, c_n])
end

function constraint_head_loss_ub_pipe_ne(wm::AbstractDirectedFlowModel,
    n::Int, a::Int, alpha::Float64, L::Float64, pipe_resistances)
    dhp = var(wm, n, :dhp, a)
    qp_ne = var(wm, n, :qp_ne, a)
    qp_ne_ub = JuMP.upper_bound.(qp_ne)
    slopes_p = pipe_resistances .* qp_ne_ub.^(alpha - 1.0)
    c_p = JuMP.@constraint(wm.model, inv(L)*dhp <= sum(slopes_p .* qp_ne))

    dhn = var(wm, n, :dhn, a)
    qn_ne = var(wm, n, :qn_ne, a)
    qn_ne_ub = JuMP.upper_bound.(qn_ne)
    slopes_n = pipe_resistances .* qn_ne_ub.^(alpha - 1.0)
    c_n = JuMP.@constraint(wm.model, inv(L)*dhn <= sum(slopes_n .* qn_ne))

    append!(con(wm, n, :head_loss)[a], [c_p, c_n])
end

function constraint_resistance_selection_ne(wm::AbstractDirectedFlowModel, n::Int, a::Int, pipe_resistances)
    c = JuMP.@constraint(wm.model, sum(var(wm, n, :x_res, a)) == 1.0)
    append!(con(wm, n, :head_loss)[a], [c])

    for (r_id, r) in enumerate(pipe_resistances)
        x_res = var(wm, n, :x_res, a)[r_id]

        qp_ne = var(wm, n, :qp_ne, a)[r_id]
        qp_ub = JuMP.upper_bound(qp_ne)
        c_p = JuMP.@constraint(wm.model, qp_ne <= qp_ub * x_res)

        qn_ne = var(wm, n, :qn_ne, a)[r_id]
        qn_ub = JuMP.upper_bound(qn_ne)
        c_n = JuMP.@constraint(wm.model, qn_ne <= qn_ub * x_res)

        append!(con(wm, n, :head_loss)[a], [c_p, c_n])
    end
end

"Constraint to ensure at least one direction is set to take flow away from a source."
function constraint_source_flow(wm::AbstractIntegerDirectedFlowForms, n::Int, i::Int, links)
    # Collect the required variables.
    x_dir = var(wm, n, :x_dir)
    out_arcs = filter(a -> i == a.second["node_fr"], links)
    out = Array{JuMP.VariableRef}([x_dir[a] for a in keys(out_arcs)])
    in_arcs = filter(a -> i == a.second["node_to"], links)
    in = Array{JuMP.VariableRef}([x_dir[a] for a in keys(in_arcs)])

    # Add the source flow direction constraint.
    c = JuMP.@constraint(wm.model, sum(out) - sum(in) >= 1.0 - length(in))
    con(wm, n, :source_flow)[i] = c
end

"Constraint to ensure at least one direction is set to take flow to a junction with demand."
function constraint_sink_flow(wm::AbstractIntegerDirectedFlowForms, n::Int, i::Int, links)
    # Collect the required variables.
    x_dir = var(wm, n, :x_dir)
    out_arcs = filter(a -> i == a.second["node_fr"], links)
    out = Array{JuMP.VariableRef}([x_dir[a] for a in keys(out_arcs)])
    in_arcs = filter(a -> i == a.second["node_to"], links)
    in = Array{JuMP.VariableRef}([x_dir[a] for a in keys(in_arcs)])

    # Add the sink flow direction constraint.
    c = JuMP.@constraint(wm.model, sum(in) - sum(out) >= 1.0 - length(out))
    con(wm, n, :sink_flow)[i] = c
end