# Constraints and variables common to all formulations with undirected flows.
# In these formulations, the variable q correspond to flow between i and j.
# When q is nonnegative, flow is assumed to travel from i to j. When q is
# negative, flow is assumed to travel from j to i.

"Create flow variables for undirected flow formulations."
function variable_flow(wm::AbstractUndirectedFlowModel; nw::Int=wm.cnw, bounded::Bool=true)
    # Create undirected flow variables (i.e., q).
    bounded ? variable_flow_bounded(wm, nw=nw) :
        variable_flow_unbounded(wm, nw=nw)
end

"Create network expansion flow variables for undirected flow formulations."
function variable_flow_ne(wm::AbstractUndirectedFlowModel; nw::Int=wm.cnw, bounded::Bool=true)
    # Create undirected flow variables (i.e., q).
    bounded ? variable_flow_bounded_ne(wm, nw=nw) :
        variable_flow_unbounded_ne(wm, nw=nw)

    # Create expressions capturing the relationships among q, and q_ne.
    var(wm, nw)[:q] = JuMP.@expression(wm.model, [a in ids(wm, nw, :link_ne)],
        sum(var(wm, nw, :q_ne, a)))

    # Create resistance binary variables.
    variable_resistance(wm, nw=nw)
end

"Create bounded flow variables for undirected flow formulations."
function variable_flow_bounded(wm::AbstractUndirectedFlowModel; nw::Int=wm.cnw)
    q_lb, q_ub = calc_flow_rate_bounds(wm, nw)

    var(wm, nw)[:q] = JuMP.@variable(wm.model, [a in ids(wm, nw, :link_fixed)],
        lower_bound=minimum(q_lb[a]), upper_bound=maximum(q_ub[a]),
        base_name="q[$(nw)]",
        start=get_start(ref(wm, nw, :link_fixed), a, "q_start", 1.0e-6))

    # Check valves and pumps are unidirectional, so flow should be nonnegative.
    lb_ids = [collect(ids(wm, nw, :check_valve)); collect(ids(wm, nw, :pump))]
    JuMP.set_lower_bound.([var(wm, nw, :q, a) for a in lb_ids], 0.0)
end

"Create bounded network expansion flow variables for undirected flow formulations."
function variable_flow_bounded_ne(wm::AbstractUndirectedFlowModel; nw::Int=wm.cnw)
    alpha = ref(wm, nw, :alpha)
    q_lb, q_ub = calc_flow_rate_bounds(wm, nw)
    link_ids = ids(wm, nw, :link_ne)
    var(wm, nw)[:q_ne] = Dict{Int,Array{JuMP.VariableRef}}(a => [] for a in link_ids)

    for a in link_ids
        n_r = length(ref(wm, nw, :resistance, a)) # Number of resistances.
        var(wm, nw, :q_ne)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            lower_bound=q_lb[a][r], upper_bound=q_ub[a][r], base_name="q_ne[$(nw)][$(a)]",
            start=get_start(ref(wm, nw, :link_ne), a, r, "q_ne_start", 1.0e-6))
    end
end

"Create unbounded flow variables for undirected flow formulations."
function variable_flow_unbounded(wm::AbstractUndirectedFlowModel; nw::Int=wm.cnw)
    var(wm, nw)[:q] = JuMP.@variable(wm.model, [a in ids(wm, nw, :link_fixed)],
        base_name="q[$(nw)]",
        start=get_start(ref(wm, nw, :link_fixed), a, "q_start", 1.0e-6))

    # Check valves and pumps are unidirectional, so flow should be nonnegative.
    lb_ids = [collect(ids(wm, nw, :check_valve)); collect(ids(wm, nw, :pump))]
    JuMP.set_lower_bound.([var(wm, nw, :q, a) for a in lb_ids], 0.0)
end

"Create unbounded network expansion flow variables for undirected flow formulations."
function variable_flow_unbounded_ne(wm::AbstractUndirectedFlowModel; nw::Int=wm.cnw)
    alpha = ref(wm, nw, :alpha)
    link_ids = ids(wm, nw, :link_ne)
    var(wm, nw)[:q_ne] = Dict{Int,Array{JuMP.VariableRef}}(a => [] for a in link_ids)

    for a in link_ids
        n_r = length(ref(wm, nw, :resistance, a)) # Number of resistances.
        var(wm, nw, :q_ne)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            base_name="q_ne[$(nw)][$(a)]",
            start=get_start(ref(wm, nw, :link_ne), a, r, "q_ne_start", 1.0e-6))
    end
end

"Constraint flow variables, based on resistance selection, in undirected flow formulations."
function constraint_resistance_selection_ne(wm::AbstractUndirectedFlowModel, n::Int, a::Int, pipe_resistances)
    c = JuMP.@constraint(wm.model, sum(var(wm, n, :x_res, a)) == 1.0)
    append!(con(wm, n, :head_loss)[a], [c])

    for r in 1:length(pipe_resistances)
        x_res = var(wm, n, :x_res, a)[r]
        q_ne = var(wm, n, :q_ne, a)[r]

        q_ne_lb = JuMP.has_lower_bound(q_ne) ? JuMP.lower_bound(q_ne) : -10.0
        c_lb = JuMP.@constraint(wm.model, q_ne >= q_ne_lb * x_res)
        q_ne_ub = JuMP.has_upper_bound(q_ne) ? JuMP.upper_bound(q_ne) : 10.0
        c_ub = JuMP.@constraint(wm.model, q_ne <= q_ne_ub * x_res)

        append!(con(wm, n, :head_loss)[a], [c_lb, c_ub])
    end
end

function constraint_check_valve(wm::AbstractUndirectedFlowModel, n::Int, a::Int, node_fr::Int, node_to::Int)
    q = var(wm, n, :q, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)
    x_cv = var(wm, n, :x_cv, a)
    q_ub = JuMP.has_upper_bound(q) ? JuMP.upper_bound(q) : 10.0

    # If the check valve is open, flow must be appreciably nonnegative.
    c_1 = JuMP.@constraint(wm.model, q <= q_ub * x_cv)
    c_2 = JuMP.@constraint(wm.model, q >= 1.0e-6 * x_cv)

    # TODO: These constraints seem to result in infeasibility in multiperiod Richmond case.
    dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    dh_ub = JuMP.upper_bound(h_i) - JuMP.lower_bound(h_j)
    c_3 = JuMP.@constraint(wm.model, h_i - h_j >= (1.0 - x_cv) * dh_lb)
    c_4 = JuMP.@constraint(wm.model, h_i - h_j <= x_cv * dh_ub)

    append!(con(wm, n, :check_valve)[a], [c_1, c_2, c_3, c_4])
end

function constraint_flow_direction_selection(wm::AbstractUndirectedFlowModel, n::Int, a::Int) end
function constraint_flow_direction_selection_ne(wm::AbstractUndirectedFlowModel, n::Int, a::Int, pipe_resistances) end
function constraint_head_difference(wm::AbstractUndirectedFlowModel, n::Int, a::Int, node_fr, node_to, head_fr, head_to)  end
function constraint_head_loss_ub_pipe_ne(wm::AbstractUndirectedFlowModel, n::Int, a::Int, alpha, len, pipe_resistances) end
function constraint_head_loss_ub_pipe(wm::AbstractUndirectedFlowModel, n::Int, a::Int, alpha, len, r_max) end
