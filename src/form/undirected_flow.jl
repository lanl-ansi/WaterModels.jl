# Constraints and variables common to all formulations with undirected flows.
# In these formulations, the variable q correspond to flow between i and j.
# When q is nonnegative, flow is assumed to travel from i to j. When q is
# negative, flow is assumed to travel from j to i.

"Create flow-related variables common to all directed flow models for edge-type components."
function variable_flow(wm::AbstractUndirectedModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    for name in ["pipe", "pump", "regulator", "short_pipe", "valve"]
        # Create undirected flow variables for each component.
        variable_component_flow(wm, name; nw=nw, bounded=bounded, report=report)
    end

    # Create flow-related variables for design components.
    variable_flow_des_common(wm, nw=nw, bounded=bounded, report=report)
end


"Create flow variables that are common to all directed flow models for a component."
function variable_component_flow(
    wm::AbstractUndirectedModel, component_name::String; nw::Int=wm.cnw,
    bounded::Bool=true, report::Bool=true)
    # Store the corresponding component symbol.
    comp_sym = Symbol(component_name)

    # Initialize the variables. (The default start value of _q_eps is crucial.)
    q = var(wm, nw)[Symbol("q_" * component_name)] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, comp_sym)], base_name="$(nw)_q_$(component_name)",
        start=comp_start_value(ref(wm, nw, comp_sym, a), "q_start", _q_eps))

    if bounded # If the variables are bounded, apply the bounds.
        q_lb, q_ub = calc_flow_bounds(wm, nw)

        for (a, comp) in ref(wm, nw, comp_sym)
            JuMP.set_lower_bound(q[a], minimum(q_lb[component_name][a]))
            JuMP.set_upper_bound(q[a], maximum(q_ub[component_name][a]))
        end
    end

    # Initialize the solution reporting data structures.
    report && sol_component_value(wm, nw, comp_sym, :q, ids(wm, nw, comp_sym), q)
end


"Create common network design flow variables for undirected flow formulations."
function variable_flow_des_common(wm::AbstractUndirectedModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    # Create dictionary for undirected design flow variables (i.e., q_des_pipe).
    q_des_pipe = var(wm, nw)[:q_des_pipe] = Dict{Int,Array{JuMP.VariableRef}}()

    # Initialize the variables. (The default start value of _q_eps is crucial.)
    for a in ids(wm, nw, :des_pipe)
        var(wm, nw, :q_des_pipe)[a] = JuMP.@variable(wm.model,
            [r in 1:length(ref(wm, nw, :resistance, a))], base_name="$(nw)_q_des_pipe",
            start=comp_start_value(ref(wm, nw, :des_pipe, a), "q_des_pipe_start", r, _q_eps))
    end

    if bounded # If the variables are bounded, apply the bounds.
        q_lb, q_ub = calc_flow_bounds(wm, nw)

        for a in ids(wm, nw, :des_pipe)
            for r in 1:length(ref(wm, nw, :resistance, a))
                JuMP.set_lower_bound(q_des_pipe[a][r], q_lb["des_pipe"][a][r])
                JuMP.set_upper_bound(q_des_pipe[a][r], q_ub["des_pipe"][a][r])
            end
        end
    end

    # Create expressions capturing the relationships among q and q_des_pipe.
    q = var(wm, nw)[:q_des_pipe_sum] = JuMP.@expression(
        wm.model, [a in ids(wm, nw, :des_pipe)], sum(var(wm, nw, :q_des_pipe, a)))

    # Initialize the solution reporting data structures.
    report && sol_component_value(wm, nw, :des_pipe, :q, ids(wm, nw, :des_pipe), q)

    # Create resistance binary variables.
    variable_resistance(wm, nw=nw)
end


function constraint_pipe_flow(wm::AbstractUndirectedModel, n::Int, a::Int)
    # By default, there are no constraints, here.
end


function constraint_pipe_head(wm::AbstractUndirectedModel, n::Int, a::Int, node_fr::Int, node_to::Int)
    # Get head variables for from and to nodes.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # For pipes, the differences must satisfy lower and upper bounds.
    dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    c_1 = JuMP.@constraint(wm.model, h_i - h_j >= dh_lb)
    dh_ub = JuMP.upper_bound(h_i) - JuMP.lower_bound(h_j)
    c_2 = JuMP.@constraint(wm.model, h_i - h_j <= dh_ub)

    # Append the constraint array.
    append!(con(wm, n, :pipe_head, a), [c_1, c_2])
end


"Constrain flow variables, based on design selections, in undirected flow formulations."
function constraint_resistance_selection_des(wm::AbstractUndirectedModel, n::Int, a::Int, pipe_resistances)
    c = JuMP.@constraint(wm.model, sum(var(wm, n, :x_res, a)) == 1.0)
    append!(con(wm, n, :head_loss)[a], [c])

    for r in 1:length(pipe_resistances)
        q_des_pipe = var(wm, n, :q_des_pipe, a)[r]
        x_res = var(wm, n, :x_res, a)[r]

        q_des_pipe_lb = JuMP.lower_bound(q_des_pipe)
        c_lb = JuMP.@constraint(wm.model, q_des_pipe >= q_des_pipe_lb * x_res)

        q_des_pipe_ub = JuMP.upper_bound(q_des_pipe)
        c_ub = JuMP.@constraint(wm.model, q_des_pipe <= q_des_pipe_ub * x_res)

        append!(con(wm, n, :head_loss)[a], [c_lb, c_ub])
    end
end




function constraint_pump_common(wm::AbstractUndirectedModel, n::Int, a::Int, node_fr::Int, node_to::Int, pc::Array{Float64})
    # Gather common variables.
    q, g, z = var(wm, n, :q_pump, a), var(wm, n, :g, a), var(wm, n, :z_pump, a)
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # If the pump is off, the flow along the pump must be zero.
    c_1 = JuMP.@constraint(wm.model, q <= JuMP.upper_bound(q) * z)
    c_2 = JuMP.@constraint(wm.model, q >= _q_eps * z)

    # If the pump is off, decouple the head difference relationship.
    dhn_lb = min(0.0, JuMP.lower_bound(h_j) - JuMP.upper_bound(h_i))
    c_3 = JuMP.@constraint(wm.model, h_j - h_i - g >= dhn_lb * (1.0 - z))
    dhn_ub = max(0.0, JuMP.upper_bound(h_j) - JuMP.lower_bound(h_i))
    c_4 = JuMP.@constraint(wm.model, h_j - h_i - g <= dhn_ub * (1.0 - z))

    # Append the constraint array.
    append!(con(wm, n, :pump, a), [c_1, c_2, c_3, c_4])
end

function constraint_pipe_common(wm::AbstractUndirectedModel, n::Int, a::Int, node_fr::Int, node_to::Int, alpha::Float64, L::Float64, r::Float64)
    # For undirected formulations, there are no constraints, here.
end


function constraint_on_off_regulator_flow(wm::AbstractUndirectedModel, n::Int, a::Int, q_min_active::Float64)
    # Get flow and regulator status variables.
    q, z = var(wm, n, :q_regulator, a), var(wm, n, :z_regulator, a)

    # If the regulator is closed, flow must be zero.
    q_lb, q_ub = max(JuMP.lower_bound(q), q_min_active), JuMP.upper_bound(q)
    c_1 = JuMP.@constraint(wm.model, q >= q_lb * z)
    c_2 = JuMP.@constraint(wm.model, q <= q_ub * z)

    # Append the constraint array.
    append!(con(wm, n, :on_off_regulator_flow, a), [c_1, c_2])
end


function constraint_on_off_regulator_head(wm::AbstractUndirectedModel, n::Int, a::Int, node_fr::Int, node_to::Int, head_setting::Float64)
    # Get regulator status variable.
    z = var(wm, n, :z_regulator, a)

    # Get head variables for from and to nodes (i.e., `i` and `j`).
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # When the pressure reducing valve is open, the head at node j is predefined.
    h_lb, h_ub = JuMP.lower_bound(h_j), JuMP.upper_bound(h_j)
    c_1 = JuMP.@constraint(wm.model, h_j >= (1.0 - z) * h_lb + z * head_setting)
    c_2 = JuMP.@constraint(wm.model, h_j <= (1.0 - z) * h_ub + z * head_setting)

    # When the pressure reducing valve is open, the head loss is nonnegative.
    dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    c_3 = JuMP.@constraint(wm.model, h_i - h_j >= dh_lb * (1.0 - z))
    dh_ub = JuMP.upper_bound(h_i) - JuMP.lower_bound(h_j)
    c_4 = JuMP.@constraint(wm.model, h_i - h_j <= dh_ub * z)

    # Append the constraint array.
    append!(con(wm, n, :on_off_regulator_head, a), [c_1, c_2, c_3, c_4])
end


function constraint_on_off_valve_flow(wm::AbstractUndirectedModel, n::Int, a::Int)
    # Get flow and valve status variables.
    q, z = var(wm, n, :q_valve, a), var(wm, n, :z_valve, a)

    # If the valve is closed, flow must be zero.
    q_lb, q_ub = JuMP.lower_bound(q), JuMP.upper_bound(q)
    c_1 = JuMP.@constraint(wm.model, q >= q_lb * z)
    c_2 = JuMP.@constraint(wm.model, q <= q_ub * z)

    # Append the constraint array.
    append!(con(wm, n, :on_off_valve_flow, a), [c_1, c_2])
end


function constraint_on_off_valve_head(wm::AbstractUndirectedModel, n::Int, a::Int, node_fr::Int, node_to::Int)
    # Get flow and valve status variables.
    q, z = var(wm, n, :q_valve, a), var(wm, n, :z_valve, a)

    # Get head variables for from and to nodes.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # When the valve is open, negative head loss is not possible.
    dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    c_1 = JuMP.@constraint(wm.model, h_i - h_j >= (1.0 - z) * dh_lb)

    # When the valve is closed, positive head loss is not possible.
    dh_ub = JuMP.upper_bound(h_i) - JuMP.lower_bound(h_j)
    c_2 = JuMP.@constraint(wm.model, h_i - h_j <= (1.0 - z) * dh_ub)

    # Append the constraint array.
    append!(con(wm, n, :on_off_valve_head, a), [c_1, c_2])
end


function constraint_on_off_pump_flow(wm::AbstractUndirectedModel, n::Int, a::Int, q_min_active::Float64)
    # Get pump status variable.
    q, z = var(wm, n, :q_pump, a), var(wm, n, :z_pump, a)

    # If the pump is inactive, flow must be zero.
    q_ub = JuMP.upper_bound(q)
    c_1 = JuMP.@constraint(wm.model, q >= q_min_active * z)
    c_2 = JuMP.@constraint(wm.model, q <= q_ub * z)

    # Append the constraint array.
    append!(con(wm, n, :on_off_pump_flow, a), [c_1, c_2])
end


function constraint_on_off_pump_head(wm::AbstractUndirectedModel, n::Int, a::Int, node_fr::Int, node_to::Int)
    # Get head variables for from and to nodes.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # Get pump status variable.
    g, z = var(wm, n, :g_pump, a), var(wm, n, :z_pump, a)

    # If the pump is off, decouple the head difference relationship. If the pump is on,
    # ensure the head difference is equal to the pump's head gain (i.e., `g`).
    dhn_ub = JuMP.upper_bound(h_j) - JuMP.lower_bound(h_i)
    dhn_lb = JuMP.lower_bound(h_j) - JuMP.upper_bound(h_i)
    c_1 = JuMP.@constraint(wm.model, h_j - h_i <= g + dhn_ub * (1.0 - z))
    c_2 = JuMP.@constraint(wm.model, h_j - h_i >= g + dhn_lb * (1.0 - z))

    # Append the constraint array.
    append!(con(wm, n, :on_off_pump_head, a), [c_1, c_2])
end


function constraint_short_pipe_flow(wm::AbstractUndirectedModel, n::Int, a::Int)
    # By default, there are no constraints, here.
end


function constraint_short_pipe_head(wm::AbstractUndirectedModel, n::Int, a::Int, node_fr::Int, node_to::Int)
    # Get head variables for from and to nodes.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # For short pipes, the heads at adjacent nodes are equal.
    c = JuMP.@constraint(wm.model, h_i - h_j == 0.0)

    # Append the constraint array.
    append!(con(wm, n, :short_pipe_head, a), [c])
end


function constraint_intermediate_directionality(
    wm::AbstractUndirectedModel, n::Int, i::Int, pipe_fr::Array{Int64,1},
    pipe_to::Array{Int64,1}, des_pipe_fr::Array{Int64,1}, des_pipe_to::Array{Int64,1},
    pump_fr::Array{Int64,1}, pump_to::Array{Int64,1}, regulator_fr::Array{Int64,1},
    regulator_to::Array{Int64,1}, short_pipe_fr::Array{Int64,1},
    short_pipe_to::Array{Int64,1}, valve_fr::Array{Int64,1}, valve_to::Array{Int64,1})
    # For undirected formulations, there are no constraints, here.
end


function constraint_sink_directionality(
    wm::AbstractUndirectedModel, n::Int, i::Int, pipe_fr::Array{Int64,1},
    pipe_to::Array{Int64,1}, des_pipe_fr::Array{Int64,1}, des_pipe_to::Array{Int64,1},
    pump_fr::Array{Int64,1}, pump_to::Array{Int64,1}, regulator_fr::Array{Int64,1},
    regulator_to::Array{Int64,1}, short_pipe_fr::Array{Int64,1},
    short_pipe_to::Array{Int64,1}, valve_fr::Array{Int64,1}, valve_to::Array{Int64,1})
    # For undirected formulations, there are no constraints, here.
end


function constraint_source_directionality(
    wm::AbstractUndirectedModel, n::Int, i::Int, pipe_fr::Array{Int64,1},
    pipe_to::Array{Int64,1}, des_pipe_fr::Array{Int64,1}, des_pipe_to::Array{Int64,1},
    pump_fr::Array{Int64,1}, pump_to::Array{Int64,1}, regulator_fr::Array{Int64,1},
    regulator_to::Array{Int64,1}, short_pipe_fr::Array{Int64,1},
    short_pipe_to::Array{Int64,1}, valve_fr::Array{Int64,1}, valve_to::Array{Int64,1})
    # For undirected formulations, there are no constraints, here.
end


function constraint_flow_direction_selection_des(wm::AbstractUndirectedModel, n::Int, a::Int, pipe_resistances) end
function constraint_pipe_head_loss_ub_des(wm::AbstractUndirectedModel, n::Int, a::Int, alpha, len, pipe_resistances) end
function constraint_pump_head_gain_lb(wm::AbstractUndirectedModel, n::Int, a::Int, node_fr::Int, node_to::Int, pc::Array{Float64}) end
