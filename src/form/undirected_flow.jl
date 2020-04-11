# Constraints and variables common to all formulations with undirected flows.
# In these formulations, the variable q correspond to flow between i and j.
# When q is nonnegative, flow is assumed to travel from i to j. When q is
# negative, flow is assumed to travel from j to i.

"Create common flow variables for undirected flow formulations."
function variable_flow_common(wm::AbstractUndirectedFlowModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
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

"Default implementation for creating flow variables for undirected flow formulations."
function variable_flow(wm::AbstractUndirectedFlowModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    variable_flow_common(wm, nw=nw, bounded=bounded, report=report)
end

"Create network design flow variables for undirected flow formulations."
function variable_flow_des(wm::AbstractUndirectedFlowModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    # Create dictionary for undirected design flow variables (i.e., q_des).
    q_des = var(wm, nw)[:q_des] = Dict{Int,Array{JuMP.VariableRef}}()

    # Initialize the variables. (The default start value of 1.0e-6 is crucial.)
    for a in ids(wm, nw, :link_des)
        var(wm, nw, :q_des)[a] = JuMP.@variable(wm.model,
            [r in 1:length(ref(wm, nw, :resistance, a))], base_name="$(nw)_q_des",
            start=comp_start_value(ref(wm, nw, :link_des, a), "q_des_start", r, 1.0e-6))
    end

    if bounded # If the variables are bounded, apply the bounds.
        q_lb, q_ub = calc_flow_bounds(wm, nw)

        for a in ids(wm, nw, :link_des)
            for r in 1:length(ref(wm, nw, :resistance, a))
                JuMP.set_lower_bound(q_des[a][r], q_lb[a][r])
                JuMP.set_upper_bound(q_des[a][r], q_ub[a][r])
            end
        end
    end

    # Create expressions capturing the relationships among q, and q_des.
    q = var(wm, nw)[:q] = JuMP.@expression(
        wm.model, [a in ids(wm, nw, :link_des)], sum(var(wm, nw, :q_des, a)))

    # Initialize the solution reporting data structures.
    report && sol_component_value(wm, nw, :link, :q, ids(wm, nw, :link_des), q)

    # Create resistance binary variables.
    variable_resistance(wm, nw=nw)
end

"Constrain flow variables, based on design selections, in undirected flow formulations."
function constraint_resistance_selection_des(wm::AbstractUndirectedFlowModel, n::Int, a::Int, pipe_resistances)
    c = JuMP.@constraint(wm.model, sum(var(wm, n, :x_res, a)) == 1.0)
    append!(con(wm, n, :head_loss)[a], [c])

    for r in 1:length(pipe_resistances)
        q_des = var(wm, n, :q_des, a)[r]
        x_res = var(wm, n, :x_res, a)[r]

        q_des_lb = JuMP.lower_bound(q_des)
        c_lb = JuMP.@constraint(wm.model, q_des >= q_des_lb * x_res)

        q_des_ub = JuMP.upper_bound(q_des)
        c_ub = JuMP.@constraint(wm.model, q_des <= q_des_ub * x_res)

        append!(con(wm, n, :head_loss)[a], [c_lb, c_ub])
    end
end

function constraint_flow_direction_selection(wm::AbstractUndirectedFlowModel, n::Int, a::Int) end
function constraint_flow_direction_selection_des(wm::AbstractUndirectedFlowModel, n::Int, a::Int, pipe_resistances) end
function constraint_head_difference(wm::AbstractUndirectedFlowModel, n::Int, a::Int, node_fr, node_to, head_fr, head_to)  end
function constraint_head_loss_ub_pipe_des(wm::AbstractUndirectedFlowModel, n::Int, a::Int, alpha, len, pipe_resistances) end
function constraint_head_loss_ub_pipe(wm::AbstractUndirectedFlowModel, n::Int, a::Int, alpha, len, r_max) end
function constraint_energy_conservation(wm::AbstractUndirectedFlowModel, n::Int, r, L, alpha) end
