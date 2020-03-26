# Constraints and variables common to all formulations with undirected flows.
# In these formulations, the variable q correspond to flow between i and j.
# When q is nonnegative, flow is assumed to travel from i to j. When q is
# negative, flow is assumed to travel from j to i.

"Create flow variables for undirected flow formulations."
function variable_flow(wm::AbstractUndirectedFlowModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    # Initialize the variables. (The default start value of 1.0e-6 is crucial.)
    q = var(wm, nw)[:q] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :link_fixed)], base_name="$(nw)_q",
        start=comp_start_value(ref(wm, nw, :link_fixed, a), "q_start", 1.0e-6))

    if bounded # If the variables are bounded, apply the bounds.
        q_lb, q_ub = calc_flow_bounds(wm, nw)

        for (a, link) in ref(wm, nw, :link_fixed)
            JuMP.set_lower_bound(q[a], minimum(q_lb[a]))
            JuMP.set_upper_bound(q[a], maximum(q_ub[a]))
        end
    end

    # Initialize the solution reporting data structures.
    report && sol_component_value(wm, nw, :link, :q, ids(wm, nw, :link_fixed), q)
end

"Create network expansion flow variables for undirected flow formulations."
function variable_flow_ne(wm::AbstractUndirectedFlowModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    # Create dictionary for undirected design flow variables (i.e., q_ne).
    q_ne = var(wm, nw)[:q_ne] = Dict{Int,Array{JuMP.VariableRef}}()

    # Initialize the variables. (The default start value of 1.0e-6 is crucial.)
    for a in ids(wm, nw, :link_ne)
        var(wm, nw, :q_ne)[a] = JuMP.@variable(wm.model,
            [r in 1:length(ref(wm, nw, :resistance, a))], base_name="$(nw)_q_ne",
            start=comp_start_value(ref(wm, nw, :link_ne, a), "q_ne_start", r, 1.0e-6))
    end

    if bounded # If the variables are bounded, apply the bounds.
        q_lb, q_ub = calc_flow_bounds(wm, nw)

        for a in ids(wm, nw, :link_ne)
            for r in 1:length(ref(wm, nw, :resistance, a))
                JuMP.set_lower_bound(q_ne[a][r], q_lb[a][r])
                JuMP.set_upper_bound(q_ne[a][r], q_ub[a][r])
            end
        end
    end

    # Create expressions capturing the relationships among q, and q_ne.
    q = var(wm, nw)[:q] = JuMP.@expression(
        wm.model, [a in ids(wm, nw, :link_ne)], sum(var(wm, nw, :q_ne, a)))

    # Initialize the solution reporting data structures.
    report && sol_component_value(wm, nw, :link, :q, ids(wm, nw, :link_ne), q)

    # Create resistance binary variables.
    variable_resistance(wm, nw=nw)
end

"Constrain flow variables, based on design selections, in undirected flow formulations."
function constraint_resistance_selection_ne(wm::AbstractUndirectedFlowModel, n::Int, a::Int, pipe_resistances)
    c = JuMP.@constraint(wm.model, sum(var(wm, n, :x_res, a)) == 1.0)
    append!(con(wm, n, :head_loss)[a], [c])

    for r in 1:length(pipe_resistances)
        q_ne = var(wm, n, :q_ne, a)[r]
        x_res = var(wm, n, :x_res, a)[r]

        q_ne_lb = JuMP.lower_bound(q_ne)
        c_lb = JuMP.@constraint(wm.model, q_ne >= q_ne_lb * x_res)

        q_ne_ub = JuMP.upper_bound(q_ne)
        c_ub = JuMP.@constraint(wm.model, q_ne <= q_ne_ub * x_res)

        append!(con(wm, n, :head_loss)[a], [c_lb, c_ub])
    end
end

function constraint_check_valve(wm::AbstractUndirectedFlowModel, n::Int, a::Int, node_fr::Int, node_to::Int)
    # Get common variables and data.
    q = var(wm, n, :q, a)
    x_cv = var(wm, n, :x_cv, a)
    h_i, h_j = [var(wm, n, :h, node_fr), var(wm, n, :h, node_to)]
    dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    dh_ub = JuMP.upper_bound(h_i) - JuMP.lower_bound(h_j)

    # If the check valve is open, flow must be appreciably nonnegative.
    c_1 = JuMP.@constraint(wm.model, q <= JuMP.upper_bound(q) * x_cv)
    c_2 = JuMP.@constraint(wm.model, q >= 6.31465679e-6 * x_cv)

    # If the check valve is open, the head difference must be nonnegative.
    c_3 = JuMP.@constraint(wm.model, h_i - h_j >= (1.0 - x_cv) * dh_lb)
    c_4 = JuMP.@constraint(wm.model, h_i - h_j <= x_cv * dh_ub)

    # Append the above to the constraint dictionary.
    append!(con(wm, n, :check_valve)[a], [c_1, c_2, c_3, c_4])
end

function constraint_flow_direction_selection(wm::AbstractUndirectedFlowModel, n::Int, a::Int) end
function constraint_flow_direction_selection_ne(wm::AbstractUndirectedFlowModel, n::Int, a::Int, pipe_resistances) end
function constraint_head_difference(wm::AbstractUndirectedFlowModel, n::Int, a::Int, node_fr, node_to, head_fr, head_to)  end
function constraint_head_loss_ub_pipe_ne(wm::AbstractUndirectedFlowModel, n::Int, a::Int, alpha, len, pipe_resistances) end
function constraint_head_loss_ub_pipe(wm::AbstractUndirectedFlowModel, n::Int, a::Int, alpha, len, r_max) end
function constraint_energy_conservation(wm::AbstractUndirectedFlowModel, n::Int, r, L, alpha) end
