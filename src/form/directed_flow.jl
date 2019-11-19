# Constraints and variables common to all formulations with directed flows. In
# these formulations, the variables qp correspond to flow from i to j, and the
# variables qn correspond to flow from j to i. That is, when qp is nonzero, qn
# should be zero, and when qn is nonzero, qp should be zero.

"Create head-related variables for formulations that use binary direction
variables. These head variables include total hydraulic head variables as well
as positive (i to j) and negative (j to i) head difference variables."
function variable_head(wm::AbstractIntegerDirectedFlowForms; nw::Int=wm.cnw, bounded::Bool=true)
    # Create variables for the total hydraulic head.
    bounded ? variable_head_bounded(wm, nw=nw) :
        variable_head_unbounded(wm, nw=nw)

    # Create variables for positive and negative hydraulic head differences.
    bounded ? variable_head_difference_bounded(wm, nw=nw) :
        variable_head_difference_unbounded(wm, nw=nw)
end

"Create bounded head difference variables for formulations with binary
direction variables. These include positive (i to j) and negative (j to i) head
difference variables, which coincide with the direction of flow."
function variable_head_difference_bounded(wm::AbstractWaterModel; nw::Int=wm.cnw)
    # Get the bounds for the head difference variables.
    dh_lb, dh_ub = calc_head_difference_bounds(wm, nw)

    # Initialize variables associated with positive head differences.
    var(wm, nw)[:dhp] = JuMP.@variable(wm.model, [a in ids(wm, nw, :link)],
        lower_bound=0.0, upper_bound=max(0.0, dh_ub[a]), base_name="dhp[$(nw)]",
        start=get_start(ref(wm, nw, :link), a, "dhp_start", 0.0))

    # Initialize variables associated with negative head differences.
    var(wm, nw)[:dhn] = JuMP.@variable(wm.model, [a in ids(wm, nw, :link)],
        lower_bound=0.0, upper_bound=max(0.0, -dh_lb[a]), base_name="dhn[$(nw)]",
        start=get_start(ref(wm, nw, :link), a, "dhn_start", 0.0))
end

"Create unbounded head difference variables for formulations with binary
direction variables. These include positive (i to j) and negative (j to i) head
difference variables, which coincide with the direction of flow."
function variable_head_difference_unbounded(wm::AbstractWaterModel; nw::Int=wm.cnw)
    # Note that the variables created below are still (very highly) bounded.
    # Initialize variables associated with positive head differences.
    var(wm, nw)[:dhp] = JuMP.@variable(wm.model, [a in ids(wm, nw, :link)],
        lower_bound=0.0, base_name="dhp[$(nw)]",
        start=get_start(ref(wm, nw, :link), a, "dhp_start", 0.0))

    # Initialize variables associated with negative head differences.
    var(wm, nw)[:dhn] = JuMP.@variable(wm.model, [a in ids(wm, nw, :link)],
        lower_bound=0.0, base_name="dhn[$(nw)]",
        start=get_start(ref(wm, nw, :link), a, "dhn_start", 0.0))
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
function variable_flow_common(wm::AbstractDirectedFlowModel; nw::Int=wm.cnw, bounded::Bool=true)
    # Create directed qp and qn variables.
    bounded ? variable_flow_bounded(wm, nw=nw) :
        variable_flow_unbounded(wm, nw=nw)

    # Create expressions capturing the relationships among q, qp, and qn.
    var(wm, nw)[:q] = JuMP.@expression(wm.model,
        [a in ids(wm, nw, :link_fixed)],
        var(wm, nw, :qp, a) - var(wm, nw, :qn, a))
end

"Create bounded flow variables for directed flow models."
function variable_flow_bounded(wm::AbstractDirectedFlowModel; nw::Int=wm.cnw)
    alpha = ref(wm, nw, :alpha)
    qn_ub, qp_ub = calc_directed_flow_upper_bounds(wm, alpha, nw)

    var(wm, nw)[:qp] = JuMP.@variable(wm.model, [a in ids(wm, nw, :link_fixed)],
        lower_bound=0.0, upper_bound=maximum(qp_ub[a]), base_name="qp[$(nw)]",
        start=get_start(ref(wm, nw, :link_fixed), a, "qp_start", 1.0e-6))

    var(wm, nw)[:qn] = JuMP.@variable(wm.model, [a in ids(wm, nw, :link_fixed)],
        lower_bound=0.0, upper_bound=maximum(qn_ub[a]), base_name="qn[$(nw)]",
        start=get_start(ref(wm, nw, :link_fixed), a, "qn_start", 1.0e-6))
end

"Create unbounded flow variables for directed flow models."
function variable_flow_unbounded(wm::AbstractDirectedFlowModel; nw::Int=wm.cnw)
    # Note that the variables created below are still bounded below by zero.
    var(wm, nw)[:qp] = JuMP.@variable(wm.model, [a in ids(wm, nw, :link_fixed)],
        lower_bound=0.0, base_name="qp[$(nw)]",
        start=get_start(ref(wm, nw, :link_fixed), a, "qp_start", 1.0e-6))

    var(wm, nw)[:qn] = JuMP.@variable(wm.model, [a in ids(wm, nw, :link_fixed)],
        lower_bound=0.0, base_name="qn[$(nw)]",
        start=get_start(ref(wm, nw, :link_fixed), a, "qn_start", 1.0e-6))
end

"Initialize variables associated with flow direction. If this variable is equal
to one, the flow direction is from i to j. If it is equal to zero, the flow
direction is from j to i."
function variable_flow_direction(wm::AbstractDirectedFlowModel; nw::Int=wm.cnw)
    #ids_bd = keys(filter(!has_check_valve, ref(wm, nw, :link)))
    var(wm, nw)[:x_dir] = JuMP.@variable(wm.model, [a in ids(wm, nw, :link)],
        base_name="x_dir[$(nw)]", binary=true, lower_bound=0, upper_bound=1,
        start=get_start(ref(wm, nw, :link), a, "x_dir_start", 1))
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
    qp = var(wm, n, :qp, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)
    x_cv = var(wm, n, :x_cv, a)

    # If the check valve is open, flow must be appreciably nonnegative.
    qp_ub = JuMP.has_upper_bound(qp) ? JuMP.upper_bound(qp) : 10.0
    c_1 = JuMP.@constraint(wm.model, qp <= qp_ub * x_cv)
    c_2 = JuMP.@constraint(wm.model, qp >= 1.0e-6 * x_cv)

    # TODO: These constraints seem to result in infeasibility in multiperiod Richmond case.
    dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    dh_ub = JuMP.upper_bound(h_i) - JuMP.lower_bound(h_j)
    c_3 = JuMP.@constraint(wm.model, h_i - h_j >= (1.0 - x_cv) * dh_lb)
    c_4 = JuMP.@constraint(wm.model, h_i - h_j <= x_cv * dh_ub)

    append!(con(wm, n, :check_valve)[a], [c_1, c_2, c_3, c_4])
end

"Create constraints for directed flow variables, based on flow direction
variables."
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
        cn = JuMP.@constraint(wm.model, qn_ne <= qn_ne_ub * (1.0 - x_dir))

        append!(con(wm, n, :head_loss)[a], [cp, cn])
    end
end

function constraint_head_difference(wm::AbstractDirectedFlowModel, n::Int, a::Int, node_fr::Int, node_to::Int, head_fr, head_to)
    x_dir = var(wm, n, :x_dir, a)

    dhp = var(wm, n, :dhp, a)
    dhp_ub = JuMP.has_upper_bound(dhp) ? JuMP.upper_bound(dhp) : 1.0e3
    cp = JuMP.@constraint(wm.model, dhp <= dhp_ub * x_dir)

    dhn = var(wm, n, :dhn, a)
    dhn_ub = JuMP.has_upper_bound(dhn) ? JuMP.upper_bound(dhn) : 1.0e3
    cn = JuMP.@constraint(wm.model, dhn <= dhn_ub * (1.0 - x_dir))

    h_i = head_fr == nothing ? var(wm, n, :h, node_fr) : head_fr
    h_j = head_to == nothing ? var(wm, n, :h, node_to) : head_to

    ce = JuMP.@constraint(wm.model, dhp - dhn == h_i - h_j)
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

    for (r_id, r) in enumerate(pipe_resistances)
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
