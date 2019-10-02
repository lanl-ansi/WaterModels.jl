# Constraints and variables common to all formulations with directed flows. In
# these formulations, the variables qp correspond to flow from i to j, and the
# variables qn correspond to flow from j to i. That is, when qp is nonzero, qn
# should be zero, and when qn is nonzero, qp should be zero.

"Create head variables for MICP, MILPR models."
function variable_head(wm::AbstractIntegerDirectedFlowForms, n::Int=wm.cnw; bounded::Bool=true)
    # Create variables for the total hydraulic head.
    bounded ? variable_head_bounded(wm, n) : variable_head_unbounded(wm, n)

    # Create variables for positive and negative hydraulic head differences.
    bounded ? variable_head_difference_bounded(wm, n) : variable_head_difference_unbounded(wm, n)
end

function variable_head_difference_bounded(wm::AbstractWaterModel, n::Int=wm.cnw)
    # Get the bounds for the head difference variables.
    dh_lb, dh_ub = calc_head_difference_bounds(wm, n)

    # Initialize variables associated with positive head differences.
    var(wm, n)[:dhp] = JuMP.@variable(wm.model, [a in ids(wm, n, :link)],
        lower_bound=0.0, upper_bound=abs(dh_ub[a]), base_name="dhp[$(n)]",
        start=get_start(ref(wm, n, :link), a, "dhp_start", 0.0))

    # Initialize variables associated with negative head differences.
    var(wm, n)[:dhn] = JuMP.@variable(wm.model, [a in ids(wm, n, :link)],
        lower_bound=0.0, upper_bound=abs(dh_lb[a]), base_name="dhn[$(n)]",
        start=get_start(ref(wm, n, :link), a, "dhn_start", 0.0))
end

function variable_head_difference_unbounded(wm::AbstractWaterModel, n::Int=wm.cnw)
    # Note that the variables created below are still (very highly) bounded.
    # Initialize variables associated with positive head differences.
    var(wm, n)[:dhp] = JuMP.@variable(wm.model, [a in ids(wm, n, :link)],
        lower_bound=0.0, base_name="dhp[$(n)]",
        start=get_start(ref(wm, n, :link), a, "dhp_start", 0.0))

    # Initialize variables associated with negative head differences.
    var(wm, n)[:dhn] = JuMP.@variable(wm.model, [a in ids(wm, n, :link)],
        lower_bound=0.0, base_name="dhn[$(n)]",
        start=get_start(ref(wm, n, :link), a, "dhn_start", 0.0))
end

"Create flow variables for CNLP model."
function variable_flow(wm::AbstractContinuousDirectedFlowForms, n::Int=wm.cnw; bounded::Bool=true)
    variable_flow_common(wm, n, bounded=bounded)
end

"Create flow variables for MICP, MILPR models."
function variable_flow(wm::AbstractIntegerDirectedFlowForms, n::Int=wm.cnw; bounded::Bool=true)
    variable_flow_common(wm, n, bounded=bounded)
    variable_flow_direction(wm, n)
end

"Create flow variables that common to all directed flow models."
function variable_flow_common(wm::AbstractDirectedFlowModel, n::Int=wm.cnw; bounded::Bool=true)
    # Create directed qp and qn variables.
    bounded ? variable_flow_bounded(wm, n) : variable_flow_unbounded(wm, n)

    # Create expressions capturing the relationships among q, qp, and qn.
    var(wm, n)[:q] = JuMP.@expression(wm.model, [a in ids(wm, n, :link_fixed)],
        var(wm, n, :qp, a) - var(wm, n, :qn, a))
end

function variable_flow_bounded(wm::AbstractDirectedFlowModel, n::Int=wm.cnw)
    alpha = ref(wm, n, :alpha)
    qn_ub, qp_ub = calc_directed_flow_upper_bounds(wm, alpha, n)

    var(wm, n)[:qp] = JuMP.@variable(wm.model, [a in ids(wm, n, :link_fixed)],
        lower_bound=0.0, upper_bound=maximum(qp_ub[a]), base_name="qp[$(n)]",
        start=get_start(ref(wm, n, :link_fixed), a, "qp_start", 1.0e-6))

    var(wm, n)[:qn] = JuMP.@variable(wm.model, [a in ids(wm, n, :link_fixed)],
        lower_bound=0.0, upper_bound=maximum(qn_ub[a]), base_name="qn[$(n)]",
        start=get_start(ref(wm, n, :link_fixed), a, "qn_start", 1.0e-6))
end

function variable_flow_unbounded(wm::AbstractDirectedFlowModel, n::Int=wm.cnw)
    # Note that the variables created below are still bounded below by zero.
    var(wm, n)[:qp] = JuMP.@variable(wm.model, [a in ids(wm, n, :link_fixed)],
        lower_bound=0.0, base_name="qp[$(n)]",
        start=get_start(ref(wm, n, :link_fixed), a, "qp_start", 1.0e-6))

    var(wm, n)[:qn] = JuMP.@variable(wm.model, [a in ids(wm, n, :link_fixed)],
        lower_bound=0.0, base_name="qn[$(n)]",
        start=get_start(ref(wm, n, :link_fixed), a, "qn_start", 1.0e-6))
end

"Initialize variables associated with flow direction. If this variable is equal
to one, the flow direction is from i to j. If it is equal to zero, the flow
direction is from j to i."
function variable_flow_direction(wm::AbstractDirectedFlowModel, n::Int=wm.cnw)
    ids_bd = keys(filter(!has_check_valve, ref(wm, n, :pipe)))
    var(wm, n)[:x_dir] = JuMP.@variable(wm.model, [a in ids_bd],
        base_name="x_dir[$(n)]", binary=true, lower_bound=0, upper_bound=1,
        start=get_start(ref(wm, n, :pipe), a, "x_dir_start", 1))
end

"Create network expansion flow variables for MICP, MILPR models."
function variable_flow_ne(wm::AbstractIntegerDirectedFlowForms, n::Int=wm.cnw; bounded::Bool=true)
    # Create directed qp_ne and qn_ne variables.
    bounded ? variable_flow_bounded_ne(wm, n) : variable_flow_unbounded_ne(wm, n)

    # Create expressions capturing the relationships among q, qp, and qn.
    var(wm, n)[:q] = JuMP.@expression(wm.model, [a in ids(wm, n, :link_ne)],
        sum(var(wm, n, :qp_ne, a)) - sum(var(wm, n, :qn_ne, a)))

    # Create resistance binary variables.
    variable_resistance(wm, n)
end

function variable_flow_ne(wm::AbstractContinuousDirectedFlowForms, n::Int=wm.cnw; bounded::Bool=true)
    if length(ref(wm, n, :link_ne)) > 0
        Memento.error(_LOGGER, "AbstractCNLPModel does not support network expansion formulations.")
    end
end

function variable_flow_bounded_ne(wm::AbstractDirectedFlowModel, n::Int=wm.cnw)
    alpha = ref(wm, n, :alpha)
    qn_ub, qp_ub = calc_directed_flow_upper_bounds(wm, alpha, n)
    link_ids = ids(wm, n, :link_ne)

    var(wm, n)[:qp_ne] = Dict{Int,Array{JuMP.VariableRef}}(a => [] for a in link_ids)
    var(wm, n)[:qn_ne] = Dict{Int,Array{JuMP.VariableRef}}(a => [] for a in link_ids)

    for a in link_ids
        n_r = length(ref(wm, n, :resistance, a)) # Number of resistances.

        var(wm, n, :qp_ne)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            lower_bound=0.0, upper_bound=qp_ub[a][r], base_name="qp_ne[$(n)][$(a)]",
            start=get_start(ref(wm, n, :link_ne), a, r, "qp_ne_start", 1.0e-6))

        var(wm, n, :qn_ne)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            lower_bound=0.0, upper_bound=qn_ub[a][r], base_name="qn_ne[$(n)][$(a)]",
            start=get_start(ref(wm, n, :link_ne), a, r, "qn_ne_start", 1.0e-6))
    end
end

function variable_flow_unbounded_ne(wm::AbstractDirectedFlowModel, n::Int=wm.cnw)
    var(wm, n)[:qp_ne] = Dict{Int,Array{JuMP.VariableRef}}
        (a => [] for a in ids(wm, n, :link_ne))
    var(wm, n)[:qn_ne] = Dict{Int,Array{JuMP.VariableRef}}
        (a => [] for a in ids(wm, n, :link_ne))

    for a in ids(wm, n, :link_ne)
        n_r = length(ref(wm, n, :resistance, a)) # Number of resistances.
        var(wm, n, :qp_ne)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            lower_bound=0.0, base_name="qp_ne[$(n)][$(a)]",
            start=get_start(ref(wm, n, :link_ne), a, "qp_ne_start", 1.0e-6))

        var(wm, n, :qn_ne)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            lower_bound=0.0, base_name="qn_ne[$(n)][$(a)]",
            start=get_start(ref(wm, n, :link_ne), a, "qn_ne_start", 1.0e-6))
    end
end

function constraint_flow_direction_selection(wm::AbstractDirectedFlowModel, n::Int, a::Int)
    x_dir = var(wm, n, :x_dir, a)

    qp = var(wm, n, :qp, a)
    qp_ub = JuMP.has_upper_bound(qp) ? JuMP.upper_bound(qp) : 10.0
    cp = JuMP.@constraint(wm.model, qp <= qp_ub * x_dir)

    qn = var(wm, n, :qn, a)
    qn_ub = JuMP.has_upper_bound(qn) ? JuMP.upper_bound(qn) : 10.0
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
        cn = JuMP.@constraint(wm.model, qn_ne - qn_ne_ub * (1.0 - x_dir) <= 0.0)

        append!(con(wm, n, :head_loss)[a], [cp, cn])
    end
end

function constraint_head_difference(wm::AbstractDirectedFlowModel, n::Int, a::Int, node_fr::Int, node_to::Int, head_fr, head_to)
    x_dir = var(wm, n, :x_dir, a)

    dhp = var(wm, n, :dhp, a)
    dhp_ub = JuMP.has_upper_bound(dhp) ? JuMP.upper_bound(dhp) : 1.0e6
    cp = JuMP.@constraint(wm.model, dhp <= dhp_ub * x_dir)

    dhn = var(wm, n, :dhn, a)
    dhn_ub = JuMP.has_upper_bound(dhn) ? JuMP.upper_bound(dhn) : 1.0e6
    cn = JuMP.@constraint(wm.model, dhn <= dhn_ub * (1.0 - x_dir))

    h_i = head_fr == nothing ? head_fr = var(wm, n, :h, node_fr) : head_fr
    h_j = head_to == nothing ? head_to = var(wm, n, :h, node_to) : head_to
    ce = JuMP.@constraint(wm.model, h_i - h_j == dhp - dhn)

    append!(con(wm, n, :head_loss)[a], [cp, cn, ce])
end

function constraint_head_loss_ub_pipe(wm::AbstractDirectedFlowModel, n::Int, a::Int, alpha::Float64, L::Float64, r::Float64)
    qp = var(wm, n, :qp, a)
    dhp = var(wm, n, :dhp, a)
    qp_ub = JuMP.has_upper_bound(qp) ? JuMP.upper_bound(qp) : 10.0
    rhs_p = r * qp_ub^(alpha - 1.0) * qp
    cp = JuMP.@constraint(wm.model, inv(L) * dhp <= rhs_p)

    qn = var(wm, n, :qn, a)
    dhn = var(wm, n, :dhn, a)
    qn_ub = JuMP.has_upper_bound(qn) ? JuMP.upper_bound(qn) : 10.0
    rhs_n = r * qn_ub^(alpha - 1.0) * qn
    cn = JuMP.@constraint(wm.model, inv(L) * dhn <= rhs_n)

    append!(con(wm, n, :head_loss)[a], [cp, cn])
end

function constraint_head_loss_ub_pipe_ne(wm::AbstractDirectedFlowModel, n::Int, a::Int, alpha::Float64, L::Float64, pipe_resistances)
    dhp = var(wm, n, :dhp, a)
    qp_ne = var(wm, n, :qp_ne, a)
    qp_ne_ub = JuMP.upper_bound.(qp_ne)
    slopes_p = pipe_resistances .* qp_ne_ub.^(alpha - 1.0)
    cp = JuMP.@constraint(wm.model, inv(L) * dhp <= sum(slopes_p .* qp_ne))

    dhn = var(wm, n, :dhn, a)
    qn_ne = var(wm, n, :qn_ne, a)
    qn_ne_ub = JuMP.upper_bound.(qn_ne)
    slopes_n = pipe_resistances .* qn_ne_ub.^(alpha - 1.0)
    cn = JuMP.@constraint(wm.model, inv(L) * dhn <= sum(slopes_n .* qn_ne))

    append!(con(wm, n, :head_loss)[a], [cp, cn])
end

function constraint_resistance_selection_ne(wm::AbstractDirectedFlowModel, n::Int, a::Int, pipe_resistances)
    c = JuMP.@constraint(wm.model, sum(var(wm, n, :x_res, a)) == 1.0)
    append!(con(wm, n, :head_loss)[a], [c])

    for r_id in 1:length(pipe_resistances)
        x_res = var(wm, n, :x_res, a)[r_id]

        qp_ne = var(wm, n, :qp_ne, a)[r_id]
        qp_ne_ub = JuMP.has_upper_bound(qp_ne) ? JuMP.upper_bound(qp_ne) : 10.0
        cp = JuMP.@constraint(wm.model, qp_ne <= qp_ne_ub * x_res)

        qn_ne = var(wm, n, :qn_ne, a)[r_id]
        qn_ne_ub = JuMP.has_upper_bound(qn_ne) ? JuMP.upper_bound(qn_ne) : 10.0
        cn = JuMP.@constraint(wm.model, qn_ne <= qn_ne_ub * x_res)

        append!(con(wm, n, :head_loss)[a], [cp, cn])
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
