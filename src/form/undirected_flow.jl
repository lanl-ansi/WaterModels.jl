# Constraints and variables common to all formulations with undirected flows.
# In these formulations, the variable q correspond to flow between i and j.
# When q is nonnegative, flow is assumed to travel from i to j. When q is
# negative, flow is assumed to travel from j to i.

function variable_flow(wm::AbstractUndirectedFlowModel, n::Int=wm.cnw; bounded::Bool=true)
    # Create undirected flow variables (i.e., q).
    bounded ? variable_flow_bounded(wm, n) : variable_flow_unbounded(wm, n)
end

function variable_flow_ne(wm::AbstractUndirectedFlowModel, n::Int=wm.cnw; bounded::Bool=true)
    # Create undirected flow variables (i.e., q).
    bounded ? variable_flow_bounded_ne(wm, n) : variable_flow_unbounded_ne(wm, n)

    # Create expressions capturing the relationships among q, and q_ne.
    var(wm, n)[:q] = JuMP.@expression(wm.model, [a in ids(wm, n, :link_ne)],
        sum(var(wm, n, :q_ne, a)))

    # Create resistance binary variables.
    variable_resistance(wm, n)
end

function variable_flow_bounded(wm::AbstractUndirectedFlowModel, n::Int=wm.cnw)
    q_lb, q_ub = calc_flow_rate_bounds(wm, n)

    var(wm, n)[:q] = JuMP.@variable(wm.model, [a in ids(wm, n, :link_fixed)],
        lower_bound=minimum(q_lb[a]), upper_bound=maximum(q_ub[a]),
        base_name="q[$(n)]",
        start=get_start(ref(wm, n, :link_fixed), a, "q_start", 1.0e-6))

    # Check valves and pumps are unidirectional, so flow should be nonnegative.
    lb_ids = [collect(ids(wm, n, :check_valve)); collect(ids(wm, n, :pump))]
    JuMP.set_lower_bound.([var(wm, n, :q, a) for a in lb_ids], 0.0)
end

function variable_flow_bounded_ne(wm::AbstractUndirectedFlowModel, n::Int=wm.cnw)
    alpha = ref(wm, n, :alpha)
    q_lb, q_ub = calc_flow_rate_bounds(wm, n)
    link_ids = ids(wm, n, :link_ne)
    var(wm, n)[:q_ne] = Dict{Int,Array{JuMP.VariableRef}}(a => [] for a in link_ids)

    for a in link_ids
        n_r = length(ref(wm, n, :resistance, a)) # Number of resistances.
        var(wm, n, :q_ne)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            lower_bound=q_lb[a][r], upper_bound=q_ub[a][r], base_name="q_ne[$(n)][$(a)]",
            start=get_start(ref(wm, n, :link_ne), a, r, "q_ne_start", 1.0e-6))
    end
end

function variable_flow_unbounded(wm::AbstractUndirectedFlowModel, n::Int=wm.cnw)
    var(wm, n)[:q] = JuMP.@variable(wm.model, [a in ids(wm, n, :link_fixed)],
        base_name="q[$(n)]",
        start=get_start(ref(wm, n, :link_fixed), a, "q_start", 1.0e-6))

    # Check valves and pumps are unidirectional, so flow should be nonnegative.
    lb_ids = [collect(ids(wm, n, :check_valve)); collect(ids(wm, n, :pump))]
    JuMP.set_lower_bound.([var(wm, n, :q, a) for a in lb_ids], 0.0)
end

function variable_flow_unbounded_ne(wm::AbstractUndirectedFlowModel, n::Int=wm.cnw)
    alpha = ref(wm, n, :alpha)
    link_ids = ids(wm, n, :link_ne)
    var(wm, n)[:q_ne] = Dict{Int,Array{JuMP.VariableRef}}(a => [] for a in link_ids)

    for a in link_ids
        n_r = length(ref(wm, n, :resistance, a)) # Number of resistances.
        var(wm, n, :q_ne)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            lower_bound=q_lb[a][r], base_name="q_ne[$(n)][$(a)]",
            start=get_start(ref(wm, n, :link_ne), a, r, "q_ne_start", 1.0e-6))
    end
end

function constraint_resistance_selection_ne(wm::AbstractUndirectedFlowModel, n::Int, a::Int, pipe_resistances)
    c = JuMP.@constraint(wm.model, sum(var(wm, n, :x_res, a)) == 1.0)
    append!(con(wm, n, :head_loss)[a], [c])

    for r in 1:length(pipe_resistances)
        x_res = var(wm, n, :x_res, a)[r]

        q_ne = var(wm, n, :q_ne, a)[r]
        q_ne_lb = JuMP.has_lower_bound(q_ne) ? JuMP.lower_bound(q_ne) : -1.0e2
        c_lb = JuMP.@constraint(wm.model, q_ne >= q_ne_lb * x_res)

        q_ne = var(wm, n, :q_ne, a)[r]
        q_ne_ub = JuMP.has_upper_bound(q_ne) ? JuMP.upper_bound(q_ne) : 1.0e2
        c_ub = JuMP.@constraint(wm.model, q_ne <= q_ne_ub * x_res)

        append!(con(wm, n, :head_loss)[a], [c_lb, c_ub])
    end
end

function constraint_flow_direction_selection(wm::AbstractUndirectedFlowModel, n::Int, a::Int) end
function constraint_flow_direction_selection_ne(wm::AbstractUndirectedFlowModel, n::Int, a::Int, pipe_resistances) end
function constraint_head_difference(wm::AbstractUndirectedFlowModel, n::Int, a::Int, node_fr, node_to, head_fr, head_to)  end
function constraint_head_loss_ub_pipe_ne(wm::AbstractUndirectedFlowModel, n::Int, a::Int, alpha, len, pipe_resistances) end
function constraint_head_loss_ub_pipe(wm::AbstractUndirectedFlowModel, n::Int, a::Int, alpha, len, r_max) end
