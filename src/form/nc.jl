# Define NC (nonconvex nonlinear programming) implementations of water distribution
# feasibility and optimization problem specifications. Note that here, ``nonconvex
# nonlinear'' describes the treatment of head loss and head gain constraints, which are
# nonlinear equalities. NC formulations also include nonlinearites that would be surmised in
# a full problem formulation (e.g., in the optimal water flow problem, a cubic function of
# flow is used in this formulation to define the consumption of power by active pumps).

# Constraints and variables common to all formulations with undirected flows.
# In these formulations, the variable q correspond to flow between i and j.
# When q is nonnegative, flow is assumed to travel from i to j. When q is
# negative, flow is assumed to travel from j to i.


"Create flow-related variables common to all directed flow models for edge-type components."
function variable_flow(
    wm::AbstractNCModel;
    nw::Int = nw_id_default,
    bounded::Bool = true,
    report::Bool = true,
)
    for name in _LINK_COMPONENTS
        # Create undirected flow variables for each component.
        _variable_component_flow(wm, name; nw = nw, bounded = bounded, report = report)
    end
end


"Create flow variables that are common to all directed flow models for a component."
function _variable_component_flow(
    wm::AbstractNCModel,
    component_name::String;
    nw::Int = nw_id_default,
    bounded::Bool = true,
    report::Bool = true,
)
    # Store the corresponding component symbol.
    comp_sym = Symbol(component_name)

    # Initialize the variables. (The default start value of _FLOW_MIN is crucial.)
    wm_data = get_wm_data(wm.data)
    flow_transform = _calc_flow_per_unit_transform(wm_data)
    flow_min_scaled = flow_transform(_FLOW_MIN)

    q =
        var(wm, nw)[Symbol("q_" * component_name)] = JuMP.@variable(
            wm.model,
            [a in ids(wm, nw, comp_sym)],
            base_name = "$(nw)_q_$(component_name)",
            start =
                comp_start_value(ref(wm, nw, comp_sym, a), "q_start", flow_min_scaled)
        )

    if bounded # If the variables are bounded, apply the bounds.
        for (a, comp) in ref(wm, nw, comp_sym)
            JuMP.set_lower_bound(q[a], comp["flow_min"])
            JuMP.set_upper_bound(q[a], comp["flow_max"])

            # Set start value for the head variable with possibly better data.
            q_mid = comp["flow_min"] + 0.5 * (comp["flow_max"] - comp["flow_min"])
            q_start_m = comp_start_value(comp, "q_start", q_mid)
            q_start = isapprox(q_start_m, 0.0; atol = 1.0e-6) ? flow_min_scaled : q_start_m
            JuMP.set_start_value(q[a], q_start)
        end
    end

    # Initialize the solution reporting data structures.
    report && sol_component_value(wm, nw, comp_sym, :q, ids(wm, nw, comp_sym), q)
end


function constraint_pipe_flow(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # By default, there are no constraints, here.
end


function constraint_pipe_head(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
)
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

"Adds head loss constraint for a pipe in the `NC` formulation."
function constraint_pipe_head_loss(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    exponent::Float64,
    L::Float64,
    r::Float64,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # Gather flow and head variables included in head loss constraints.
    q, h_i, h_j = var(wm, n, :q_pipe, a), var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # Add nonconvex constraint for the head loss relationship.
    c_1 = JuMP.@NLconstraint(wm.model, r * head_loss(q) <= (h_i - h_j) / L)
    c_2 = JuMP.@NLconstraint(wm.model, r * head_loss(q) >= (h_i - h_j) / L)

    # Append the :pipe_head_loss constraint array.
    append!(con(wm, n, :pipe_head_loss)[a], [c_1, c_2])
end


function constraint_on_off_des_pipe_flow(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # Get flow and design status variables.
    q, z = var(wm, n, :q_des_pipe, a), var(wm, n, :z_des_pipe, a)

    # If the valve is closed, flow must be zero.
    q_lb, q_ub = JuMP.lower_bound(q), JuMP.upper_bound(q)
    c_1 = JuMP.@constraint(wm.model, q >= q_lb * z)
    c_2 = JuMP.@constraint(wm.model, q <= q_ub * z)

    # Append the constraint array.
    append!(con(wm, n, :on_off_des_pipe_flow, a), [c_1, c_2])
end


function constraint_on_off_des_pipe_head(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
)
    # Get head difference and status variables for the design pipe.
    dh = var(wm, n, :dh_des_pipe, a)
    z = var(wm, n, :z_des_pipe, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)

    # For pipes, the differences must satisfy lower and upper bounds.
    dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    c_1 = JuMP.@constraint(wm.model, dh >= dh_lb * z)
    dh_ub = JuMP.upper_bound(h_i) - JuMP.lower_bound(h_j)
    c_2 = JuMP.@constraint(wm.model, dh <= dh_ub * z)

    # Append the :n_off_des_pipe_head constraint array.
    append!(con(wm, n, :on_off_des_pipe_head, a), [c_1, c_2])
end


"Adds head loss constraint for a design pipe in the `NC` formulation."
function constraint_on_off_des_pipe_head_loss(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    exponent::Float64,
    L::Float64,
    r::Float64,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # Get flow and head difference variables.
    q, dh = var(wm, n, :q_des_pipe, a), var(wm, n, :dh_des_pipe, a)

    # Add nonconvex constraint for the head loss relationship.
    c_1 = JuMP.@NLconstraint(wm.model, r * head_loss(q) <= dh / L)
    c_2 = JuMP.@NLconstraint(wm.model, r * head_loss(q) >= dh / L)

    # Append the :pipe_head_loss constraint array.
    append!(con(wm, n, :on_off_des_pipe_head_loss)[a], [c_1, c_2])
end


function constraint_des_pipe_flow(
    wm::AbstractNCModel,
    n::Int,
    k::Int,
    node_fr::Int,
    node_to::Int,
    des_pipes::Vector{Int},
)
    # For undirected formulations, there are no constraints, here.
end


function constraint_des_pipe_head(
    wm::AbstractNCModel,
    n::Int,
    k::Int,
    node_fr::Int,
    node_to::Int,
    des_pipes::Vector{Int},
)
    # Get head-related variables for design pipes and arcs.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)
    dh_sum = sum(var(wm, n, :dh_des_pipe, a) for a in des_pipes)

    # Add constraint equating the sum of head differences to head difference.
    c = JuMP.@constraint(wm.model, dh_sum == h_i - h_j)

    # Append the :pipe_head_loss constraint array.
    append!(con(wm, n, :des_pipe_head)[k], [c])
end


function constraint_on_off_pump_flow(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    q_min_forward::Float64,
)
    # Get pump status variable.
    q, z = var(wm, n, :q_pump, a), var(wm, n, :z_pump, a)

    # If the pump is inactive, flow must be zero.
    q_ub = JuMP.upper_bound(q)
    c_1 = JuMP.@constraint(wm.model, q >= q_min_forward * z)
    c_2 = JuMP.@constraint(wm.model, q <= q_ub * z)

    # Append the constraint array.
    append!(con(wm, n, :on_off_pump_flow, a), [c_1, c_2])
end

function constraint_on_off_pump_flow_ne(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    q_min_forward::Float64,
)
    if(n==1)
        println("Running NC pump flow ne")
    end
    # Get pump status variable.
    q, z = var(wm, n, :q_ne_pump, a), var(wm, n, :z_ne_pump, a)

    # If the pump is inactive, flow must be zero.
    q_ub = JuMP.upper_bound(q)
    c_1 = JuMP.@constraint(wm.model, q >= q_min_forward * z)
    c_2 = JuMP.@constraint(wm.model, q <= q_ub * z)

    # Append the constraint array.
    append!(con(wm, n, :on_off_pump_flow_ne, a), [c_1, c_2])
end

function constraint_on_off_pump_head(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
)
    if(n==1)
        println("Running nc pump head constraint")
    end
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

function constraint_on_off_pump_head_ne(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
)
    if(n == 1)
        println("Running NC pump head ne")
    end
    # Get head variables for from and to nodes.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # Get pump status variable.
    g, z = var(wm, n, :g_ne_pump, a), var(wm, n, :z_ne_pump, a)

    # If the pump is off, decouple the head difference relationship. If the pump is on,
    # ensure the head difference is equal to the pump's head gain (i.e., `g`).
    dhn_ub = JuMP.upper_bound(h_j) - JuMP.lower_bound(h_i)
    dhn_lb = JuMP.lower_bound(h_j) - JuMP.upper_bound(h_i)
    c_1 = JuMP.@constraint(wm.model, h_j - h_i <= g + dhn_ub * (1.0 - z))
    c_2 = JuMP.@constraint(wm.model, h_j - h_i >= g + dhn_lb * (1.0 - z))

    # Append the constraint array.
    append!(con(wm, n, :on_off_pump_head_ne, a), [c_1, c_2])
end

"Adds head gain constraints for pumps in `NC` formulations."
function constraint_on_off_pump_head_gain(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    q_min_forward::Float64,
)
    # Gather pump flow, head gain, and status variables.
    q, g, z = var(wm, n, :q_pump, a), var(wm, n, :g_pump, a), var(wm, n, :z_pump, a)

    # Add constraint equating head gain with respect to the pump curve.
    head_curve_func = _calc_head_curve_function(ref(wm, n, :pump, a), z)
    c_1 = JuMP.@constraint(wm.model, head_curve_func(q) <= g)
    c_2 = JuMP.@constraint(wm.model, head_curve_func(q) >= g)

    # Append the :on_off_pump_head_gain constraint array.
    append!(con(wm, n, :on_off_pump_head_gain)[a], [c_1, c_2])
end

"Adds head gain constraints for expansion pumps in `NC` formulations."
function constraint_on_off_pump_head_gain_ne(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    q_min_forward::Float64,
)
    if(n==1)
        println("Running NC pump head gain ne")
    end
    # Gather pump flow, head gain, and status variables.
    q, g, z = var(wm, n, :q_ne_pump, a), var(wm, n, :g_ne_pump, a), var(wm, n, :z_ne_pump, a)

    # Add constraint equating head gain with respect to the pump curve.
    head_curve_func = _calc_head_curve_function(ref(wm, n, :ne_pump, a), z)
    c_1 = JuMP.@constraint(wm.model, head_curve_func(q) <= g)
    c_2 = JuMP.@constraint(wm.model, head_curve_func(q) >= g)

    # Append the :on_off_pump_head_gain constraint array.
    append!(con(wm, n, :on_off_pump_head_gain_ne)[a], [c_1, c_2])
end

function constraint_on_off_pump_power(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    q_min_forward::Float64,
)
    # Gather pump flow, scaled power, and status variables.
    q, P, z = var(wm, n, :q_pump, a), var(wm, n, :P_pump, a), var(wm, n, :z_pump, a)

    # Add constraint equating power with respect to the power curve.
    power_qa = _calc_pump_power_quadratic_approximation(wm, n, a, z)
    c_1 = JuMP.@constraint(wm.model, power_qa(q) <= P)
    c_2 = JuMP.@constraint(wm.model, power_qa(q) >= P)

    # Append the :on_off_pump_power constraint array.
    append!(con(wm, n, :on_off_pump_power)[a], [c_1, c_2])
end

function constraint_on_off_pump_power_ne(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    q_min_forward::Float64,
)
    if(n == 1)
        println("Running NC pump power ne")
    end
    # Gather pump flow, scaled power, and status variables.
    q, P, z = var(wm, n, :q_ne_pump, a), var(wm, n, :P_ne_pump, a), var(wm, n, :z_ne_pump, a)

    # Add constraint equating power with respect to the power curve.
    power_qa = _calc_pump_power_quadratic_approximation_ne(wm, n, a, z)
    c_1 = JuMP.@constraint(wm.model, power_qa(q) <= P)
    c_2 = JuMP.@constraint(wm.model, power_qa(q) >= P)

    # Append the :on_off_pump_power constraint array.
    append!(con(wm, n, :on_off_pump_power_ne)[a], [c_1, c_2])
end


function constraint_on_off_regulator_flow(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    q_min_forward::Float64,
)
    # Get flow and regulator status variables.
    q, z = var(wm, n, :q_regulator, a), var(wm, n, :z_regulator, a)

    # If the regulator is closed, flow must be zero.
    q_lb, q_ub = max(JuMP.lower_bound(q), q_min_forward), JuMP.upper_bound(q)
    c_1 = JuMP.@constraint(wm.model, q >= q_lb * z)
    c_2 = JuMP.@constraint(wm.model, q <= q_ub * z)

    # Append the constraint array.
    append!(con(wm, n, :on_off_regulator_flow, a), [c_1, c_2])
end


function constraint_on_off_regulator_head(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    head_setting::Float64,
)
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


function constraint_on_off_valve_flow(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # Get flow and valve status variables.
    q, z = var(wm, n, :q_valve, a), var(wm, n, :z_valve, a)

    # If the valve is closed, flow must be zero.
    q_lb, q_ub = JuMP.lower_bound(q), JuMP.upper_bound(q)
    c_1 = JuMP.@constraint(wm.model, q >= q_lb * z)
    c_2 = JuMP.@constraint(wm.model, q <= q_ub * z)

    # Append the constraint array.
    append!(con(wm, n, :on_off_valve_flow, a), [c_1, c_2])
end


function constraint_on_off_valve_head(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
)
    # Get flow and valve status variables.
    q, z = var(wm, n, :q_valve, a), var(wm, n, :z_valve, a)

    # Get head variables for from and to nodes.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # When the valve is open, negative head loss is not possible.
    dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    c_1 = JuMP.@constraint(wm.model, h_i - h_j >= (1.0 - z) * dh_lb)

    # When the valve is open, positive head loss is not possible.
    dh_ub = JuMP.upper_bound(h_i) - JuMP.lower_bound(h_j)
    c_2 = JuMP.@constraint(wm.model, h_i - h_j <= (1.0 - z) * dh_ub)

    # Append the constraint array.
    append!(con(wm, n, :on_off_valve_head, a), [c_1, c_2])
end


function constraint_short_pipe_flow(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # By default, there are no constraints, here.
end


function constraint_short_pipe_flow_ne(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # Get flow and status variables for the short pipe.
    q = var(wm, n, :q_ne_short_pipe, a)
    z = var(wm, n, :z_ne_short_pipe, a)

    # Get lower and upper bounds for the expansion short pipe's flow.
    q_lb = JuMP.lower_bound(q)
    q_ub = JuMP.upper_bound(q)

    # Add constraints limiting the flow based on expansion status.
    c_1 = JuMP.@constraint(wm.model, q <= q_ub * z)
    c_2 = JuMP.@constraint(wm.model, q >= q_lb * z)

    # Append the constraint array.
    append!(con(wm, n, :short_pipe_flow_ne, a), [c_1, c_2])
end


function constraint_short_pipe_head(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
)
    # Get head variables for from and to nodes.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # For short pipes, the heads at adjacent nodes are equal.
    c = JuMP.@constraint(wm.model, h_i == h_j)

    # Append the constraint array.
    append!(con(wm, n, :short_pipe_head, a), [c])
end


function constraint_short_pipe_head_ne(
    wm::AbstractNCModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
)
    # Get expansion shot pipe status variable.
    z = var(wm, n, :z_ne_short_pipe, a)

    # Get head variables for from and to nodes.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # When the short pipe is constructed, negative head loss is not possible.
    dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    c_1 = JuMP.@constraint(wm.model, h_i - h_j >= (1.0 - z) * dh_lb)

    # When the short pipe is constructed, positive head loss is not possible.
    dh_ub = JuMP.upper_bound(h_i) - JuMP.lower_bound(h_j)
    c_2 = JuMP.@constraint(wm.model, h_i - h_j <= (1.0 - z) * dh_ub)

    # Append the constraint array.
    append!(con(wm, n, :short_pipe_head_ne, a), [c_1, c_2])
end


function constraint_intermediate_directionality(
    wm::AbstractNCModel,
    n::Int,
    i::Int,
    pipe_fr::Vector{Int},
    pipe_to::Vector{Int},
    des_pipe_fr::Vector{Int},
    des_pipe_to::Vector{Int},
    pump_fr::Vector{Int},
    pump_to::Vector{Int},
    ne_pump_fr::Vector{Int},
    ne_pump_to::Vector{Int},
    regulator_fr::Vector{Int},
    regulator_to::Vector{Int},
    short_pipe_fr::Vector{Int},
    short_pipe_to::Vector{Int},
    ne_short_pipe_fr::Vector{Int},
    ne_short_pipe_to::Vector{Int},
    valve_fr::Vector{Int},
    valve_to::Vector{Int},
)
    # For undirected formulations, there are no constraints, here.
end


function constraint_sink_directionality(
    wm::AbstractNCModel,
    n::Int,
    i::Int,
    pipe_fr::Vector{Int},
    pipe_to::Vector{Int},
    des_pipe_fr::Vector{Int},
    des_pipe_to::Vector{Int},
    pump_fr::Vector{Int},
    pump_to::Vector{Int},
    ne_pump_fr::Vector{Int},
    ne_pump_to::Vector{Int},
    regulator_fr::Vector{Int},
    regulator_to::Vector{Int},
    short_pipe_fr::Vector{Int},
    short_pipe_to::Vector{Int},
    ne_short_pipe_fr::Vector{Int},
    ne_short_pipe_to::Vector{Int},
    valve_fr::Vector{Int},
    valve_to::Vector{Int},
)
    # For undirected formulations, there are no constraints, here.
end


function constraint_source_directionality(
    wm::AbstractNCModel,
    n::Int,
    i::Int,
    pipe_fr::Vector{Int},
    pipe_to::Vector{Int},
    des_pipe_fr::Vector{Int},
    des_pipe_to::Vector{Int},
    pump_fr::Vector{Int},
    pump_to::Vector{Int},
    ne_pump_fr::Vector{Int},
    ne_pump_to::Vector{Int},
    regulator_fr::Vector{Int},
    regulator_to::Vector{Int},
    short_pipe_fr::Vector{Int},
    short_pipe_to::Vector{Int},
    ne_short_pipe_fr::Vector{Int},
    ne_short_pipe_to::Vector{Int},
    valve_fr::Vector{Int},
    valve_to::Vector{Int},
)
    # For undirected formulations, there are no constraints, here.
end
