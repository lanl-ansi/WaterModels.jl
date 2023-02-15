# Constraint templates help simplify data wrangling across multiple water network
# optimization formulations by providing an abstraction layer between the network data and
# network constraint definitions. Each constraint template's job is to extract the required
# parameters from a given network data structure and pass the data as named arguments to the
# water network optimization constraints. Here, constraint templates should always be
# defined over the type `AbstractWaterModel` and should never refer to model variables.


function _initialize_con_dict(
    wm::AbstractWaterModel,
    key::Symbol;
    nw::Int = nw_id_default,
    is_array::Bool = false,
)
    if !haskey(con(wm, nw), key)
        ref_type = is_array ? Vector{JuMP.ConstraintRef} : JuMP.ConstraintRef
        con(wm, nw)[key] = Dict{Int,ref_type}()
    end
end


function _get_from_and_to_components(
    wm::AbstractWaterModel,
    i::Int,
    component_type::Symbol;
    nw::Int = nw_id_default,
)
    fr_sym = Symbol(string(component_type) * "_fr")
    to_sym = Symbol(string(component_type) * "_to")
    return ref(wm, nw, fr_sym, i), ref(wm, nw, to_sym, i)
end


### Nodal Constraints ###
function constraint_flow_conservation(
    wm::AbstractWaterModel,
    i::Int;
    nw::Int = nw_id_default,
)
    # Collect various indices for edge-type components connected to node `i`.
    pipe_fr, pipe_to = _get_from_and_to_components(wm, i, :pipe; nw = nw)
    des_pipe_fr, des_pipe_to = _get_from_and_to_components(wm, i, :des_pipe; nw = nw)
    pump_fr, pump_to = _get_from_and_to_components(wm, i, :pump; nw = nw)
    ne_pump_fr, ne_pump_to = _get_from_and_to_components(wm, i, :ne_pump; nw = nw)
    regulator_fr, regulator_to = _get_from_and_to_components(wm, i, :regulator; nw = nw)
    short_pipe_fr, short_pipe_to = _get_from_and_to_components(wm, i, :short_pipe; nw = nw)
    ne_short_pipe_fr, ne_short_pipe_to = _get_from_and_to_components(wm, i, :ne_short_pipe; nw = nw)
    valve_fr, valve_to = _get_from_and_to_components(wm, i, :valve; nw = nw)

    # Collect various indices for node-type components connected to node `i`.
    reservoirs = ref(wm, nw, :node_reservoir, i) # Reservoirs attached to node `i`.
    tanks = ref(wm, nw, :node_tank, i) # Tanks attached to node `i`.
    demands = ref(wm, nw, :node_demand, i) # Demands attached to node `i`.

    # Sum the constant demands required at node `i`.
    nondispatch_demand_ids = ids(wm, nw, :nondispatchable_demand)
    nd_demands_at_i = filter(j -> j in nondispatch_demand_ids, demands)
    fixed_demands_at_i = ref.(Ref(wm), nw, :nondispatchable_demand, nd_demands_at_i)
    fixed_flows_at_i = [x["flow_nominal"] for x in fixed_demands_at_i]
    net_fixed_demand = length(fixed_flows_at_i) > 0 ? sum(fixed_flows_at_i) : 0.0

    # Get the indices of dispatchable demands connected to node `i`.
    dispatchable_demands = filter(j -> j in ids(wm, nw, :dispatchable_demand), demands)

    # Initialize the flow conservation constraint dictionary entry.
    _initialize_con_dict(wm, :flow_conservation, nw = nw)

    # Add the flow conservation constraint.
    constraint_flow_conservation(
        wm,
        nw,
        i,
        pipe_fr,
        pipe_to,
        des_pipe_fr,
        des_pipe_to,
        pump_fr,
        pump_to,
        ne_pump_fr,
        ne_pump_to,
        regulator_fr,
        regulator_to,
        short_pipe_fr,
        short_pipe_to,
        ne_short_pipe_fr,
        ne_short_pipe_to,
        valve_fr,
        valve_to,
        reservoirs,
        tanks,
        dispatchable_demands,
        net_fixed_demand,
    )
end


function constraint_node_directionality(
    wm::AbstractWaterModel,
    i::Int;
    nw::Int = nw_id_default,
)
    # Collect various indices for edge-type components connected to node `i`.
    pipe_fr, pipe_to = _get_from_and_to_components(wm, i, :pipe; nw = nw)
    des_pipe_fr, des_pipe_to = _get_from_and_to_components(wm, i, :des_pipe; nw = nw)
    pump_fr, pump_to = _get_from_and_to_components(wm, i, :pump; nw = nw)
    ne_pump_fr, ne_pump_to = _get_from_and_to_components(wm, i, :ne_pump; nw = nw)
    regulator_fr, regulator_to = _get_from_and_to_components(wm, i, :regulator; nw = nw)
    short_pipe_fr, short_pipe_to = _get_from_and_to_components(wm, i, :short_pipe; nw = nw)
    ne_short_pipe_fr, ne_short_pipe_to = _get_from_and_to_components(wm, i, :ne_short_pipe; nw = nw)
    valve_fr, valve_to = _get_from_and_to_components(wm, i, :valve; nw = nw)

    # Collect various indices for node-type components connected to node `i`.
    reservoirs = ref(wm, nw, :node_reservoir, i) # Reservoirs attached to node `i`.
    tanks = ref(wm, nw, :node_tank, i) # Tanks attached to node `i`.
    demands = ref(wm, nw, :node_demand, i) # Demands attached to node `i`.

    # Sum the constant demands required at node `i`.
    nondispatch_demand_ids = ids(wm, nw, :nondispatchable_demand)
    nd_demands_at_i = filter(j -> j in nondispatch_demand_ids, demands)
    fixed_demands_at_i = ref.(Ref(wm), nw, :nondispatchable_demand, nd_demands_at_i)
    fixed_flows_at_i = [x["flow_nominal"] for x in fixed_demands_at_i]
    net_fixed_demand = length(fixed_flows_at_i) > 0 ? sum(fixed_flows_at_i) : 0.0

    # Get the in degree of node `i`.
    in_length =
        length(pipe_to) +
        length(des_pipe_to) +
        length(pump_to) +
        length(ne_pump_to) +
        length(regulator_to) +
        length(short_pipe_to) +
        length(ne_short_pipe_to) +
        length(valve_to)

    # Get the out degree of node `i`.
    out_length =
        length(pipe_fr) +
        length(des_pipe_fr) +
        length(pump_fr) +
        length(ne_pump_fr) +
        length(regulator_fr) +
        length(short_pipe_fr) +
        length(ne_short_pipe_fr) +
        length(valve_fr)

    # Initialize the directionality constraint dictionary entry.
    _initialize_con_dict(wm, :node_directionality, nw = nw)

    # Get the number of nodal components attached to node `i`.
    num_components = length(demands) + length(tanks) + length(reservoirs)

    # Check if node directionality constraints should be added.
    if num_components == 0 && in_length + out_length == 2
        # Add the intermediate node directionality constraint.
        constraint_intermediate_directionality(
            wm,
            nw,
            i,
            pipe_fr,
            pipe_to,
            des_pipe_fr,
            des_pipe_to,
            pump_fr,
            pump_to,
            ne_pump_fr,
            ne_pump_to,
            regulator_fr,
            regulator_to,
            short_pipe_fr,
            short_pipe_to,
            ne_short_pipe_fr,
            ne_short_pipe_to,
            valve_fr,
            valve_to,
        )
    elseif length(reservoirs) > 0 && num_components == length(reservoirs)
        # Add the source node directionality constraint.
        constraint_source_directionality(
            wm,
            nw,
            i,
            pipe_fr,
            pipe_to,
            des_pipe_fr,
            des_pipe_to,
            pump_fr,
            pump_to,
            ne_pump_fr,
            ne_pump_to,
            regulator_fr,
            regulator_to,
            short_pipe_fr,
            short_pipe_to,
            ne_short_pipe_fr,
            ne_short_pipe_to,
            valve_fr,
            valve_to,
        )
    elseif num_components == length(fixed_demands_at_i) && net_fixed_demand < 0.0
        # Add the source node directionality constraint.
        constraint_source_directionality(
            wm,
            nw,
            i,
            pipe_fr,
            pipe_to,
            des_pipe_fr,
            des_pipe_to,
            pump_fr,
            pump_to,
            ne_pump_fr,
            ne_pump_to,
            regulator_fr,
            regulator_to,
            short_pipe_fr,
            short_pipe_to,
            ne_short_pipe_fr,
            ne_short_pipe_to,
            valve_fr,
            valve_to,
        )
    elseif num_components == length(fixed_demands_at_i) && net_fixed_demand > 0.0
        # Add the sink node directionality constraint.
        constraint_sink_directionality(
            wm,
            nw,
            i,
            pipe_fr,
            pipe_to,
            des_pipe_fr,
            des_pipe_to,
            pump_fr,
            pump_to,
            ne_pump_fr,
            ne_pump_to,
            regulator_fr,
            regulator_to,
            short_pipe_fr,
            short_pipe_to,
            ne_short_pipe_fr,
            ne_short_pipe_to,
            valve_fr,
            valve_to,
        )
    end
end


### Tank Constraints ###
function constraint_tank_volume(wm::AbstractWaterModel, i::Int; nw::Int = nw_id_default)
    if !ref(wm, nw, :tank, i, "dispatchable")
        # Only set the tank state if the tank is nondispatchable.
        constraint_tank_volume_fixed(wm, i; nw = nw)
    end
end


function constraint_tank_volume_fixed(
    wm::AbstractWaterModel,
    i::Int;
    nw::Int = nw_id_default,
)
    # Compute the cross-sectional area of the tank.
    tank = ref(wm, nw, :tank, i)
    initial_level = tank["init_level"]
    surface_area = 0.25 * pi * tank["diameter"]^2

    # Compute initial, minimum, and maximum tank volumes.
    V_initial = surface_area * initial_level
    time_step = ref(wm, nw, :time_step)
    V_min = max(tank["min_vol"], surface_area * tank["min_level"])
    V_max = surface_area * tank["max_level"]

    # Update the nodal elevation data corresponding to the tank.
    node = ref(wm, nw, :node, tank["node"])
    head = node["elevation"] + initial_level
    node["head_nominal"] = node["head_min"] = node["head_max"] = head

    # Apply the tank volume constraint at the specified time step.
    _initialize_con_dict(wm, :tank_volume; nw = nw, is_array = true)
    con(wm, nw, :tank_volume)[i] = Array{JuMP.ConstraintRef,1}([])
    constraint_tank_volume_fixed(wm, nw, i, V_initial, time_step, V_min, V_max)
end


function constraint_tank_volume(wm::AbstractWaterModel, i::Int, nw_1::Int, nw_2::Int)
    # Only apply the constraint if the tank exists in both subnetworks.
    if haskey(ref(wm, nw_1, :tank), i) && haskey(ref(wm, nw_2, :tank), i)
        # Get the tank reference within each of the subnetworks.
        tank_nw_1 = ref(wm, nw_1, :tank, i)
        tank_nw_2 = ref(wm, nw_2, :tank, i)

        # Only set the tank state if the tank is nondispatchable.
        if !tank_nw_1["dispatchable"] && !tank_nw_2["dispatchable"]
            # Apply the tank volume integration constraint between the two time steps.
            _initialize_con_dict(wm, :tank_volume; nw = nw_2, is_array = true)
            con(wm, nw_2, :tank_volume)[i] = Array{JuMP.ConstraintRef,1}([])
            constraint_tank_volume(wm, nw_1, nw_2, i, ref(wm, nw_1, :time_step))
        end
    end
end


### Pipe Constraints ###
function constraint_pipe_flow(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    pipe = ref(wm, nw, :pipe, a)
    q_max_reverse = min(get(pipe, "flow_max_reverse", 0.0), pipe["flow_max"])
    q_min_forward = max(get(pipe, "flow_min_forward", 0.0), pipe["flow_min"])

    _initialize_con_dict(wm, :pipe_flow, nw = nw, is_array = true)
    con(wm, nw, :pipe_flow)[a] = Array{JuMP.ConstraintRef}([])
    constraint_pipe_flow(wm, nw, a, q_max_reverse, q_min_forward)
end


function constraint_pipe_head(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr = ref(wm, nw, :pipe, a, "node_fr")
    node_to = ref(wm, nw, :pipe, a, "node_to")
    _initialize_con_dict(wm, :pipe_head, nw = nw, is_array = true)
    con(wm, nw, :pipe_head)[a] = Array{JuMP.ConstraintRef}([])
    constraint_pipe_head(wm, nw, a, node_fr, node_to)
end


function constraint_pipe_head_loss(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr = ref(wm, nw, :pipe, a, "node_fr")
    node_to = ref(wm, nw, :pipe, a, "node_to")
    exponent = _get_exponent_from_head_loss_form(wm.ref[:it][wm_it_sym][:head_loss])
    L = ref(wm, nw, :pipe, a, "length")

    wm_data = get_wm_data(wm.data)
    head_loss, viscosity = wm_data["head_loss"], wm_data["viscosity"]
    base_length = get(wm_data, "base_length", 1.0)
    base_mass = get(wm_data, "base_mass", 1.0)
    base_time = get(wm_data, "base_time", 1.0)

    r = _calc_pipe_resistance(
        ref(wm, nw, :pipe, a),
        head_loss,
        viscosity,
        base_length,
        base_mass,
        base_time,
    )
    q_max_reverse = min(get(ref(wm, nw, :pipe, a), "flow_max_reverse", 0.0), 0.0)
    q_min_forward = max(get(ref(wm, nw, :pipe, a), "flow_min_forward", 0.0), 0.0)

    _initialize_con_dict(wm, :pipe_head_loss, nw = nw, is_array = true)
    con(wm, nw, :pipe_head_loss)[a] = Array{JuMP.ConstraintRef}([])
    constraint_pipe_head_loss(
        wm,
        nw,
        a,
        node_fr,
        node_to,
        exponent,
        L,
        r,
        q_max_reverse,
        q_min_forward,
    )
end


### Design Pipe Constraints ###
function constraint_des_pipe_flow(
    wm::AbstractWaterModel,
    k::Int,
    node_fr::Int,
    node_to::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    des_pipes = collect(
        keys(
            filter(
                x -> x.second["node_fr"] == node_fr && x.second["node_to"] == node_to,
                ref(wm, nw, :des_pipe),
            ),
        ),
    )

    _initialize_con_dict(wm, :des_pipe_flow, nw = nw, is_array = true)
    con(wm, nw, :des_pipe_flow)[k] = Array{JuMP.ConstraintRef}([])
    constraint_des_pipe_flow(wm, nw, k, node_fr, node_to, des_pipes)
end


function constraint_des_pipe_head(
    wm::AbstractWaterModel,
    k::Int,
    node_fr::Int,
    node_to::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    des_pipes = collect(
        keys(
            filter(
                x -> x.second["node_fr"] == node_fr && x.second["node_to"] == node_to,
                ref(wm, nw, :des_pipe),
            ),
        ),
    )

    _initialize_con_dict(wm, :des_pipe_head, nw = nw, is_array = true)
    con(wm, nw, :des_pipe_head)[k] = Array{JuMP.ConstraintRef}([])
    constraint_des_pipe_head(wm, nw, k, node_fr, node_to, des_pipes)
end


function constraint_des_pipe_selection(
    wm::AbstractWaterModel,
    k::Int,
    node_fr::Int,
    node_to::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    des_pipes = collect(
        keys(
            filter(
                x -> x.second["node_fr"] == node_fr && x.second["node_to"] == node_to,
                ref(wm, nw, :des_pipe),
            ),
        ),
    )

    _initialize_con_dict(wm, :des_pipe_selection, nw = nw, is_array = true)
    con(wm, nw, :des_pipe_selection)[k] = Array{JuMP.ConstraintRef}([])
    constraint_des_pipe_selection(wm, nw, k, node_fr, node_to, des_pipes)
end


function constraint_on_off_des_pipe_flow(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    # Get the design pipe from the specified index.
    des_pipe = ref(wm, nw, :des_pipe, a)

    # Compute metadata associated with the design pipe.
    q_max_reverse = min(get(des_pipe, "flow_max_reverse", 0.0), 0.0)
    q_min_forward = max(get(des_pipe, "flow_min_forward", 0.0), 0.0)

    # Initialize :on_off_des_pipe_flow constraint dictionary.
    _initialize_con_dict(wm, :on_off_des_pipe_flow, nw = nw, is_array = true)
    con(wm, nw, :on_off_des_pipe_flow)[a] = Array{JuMP.ConstraintRef}([])

    # Apply the :on_off_des_pipe_flow constraints.
    constraint_on_off_des_pipe_flow(wm, nw, a, q_max_reverse, q_min_forward)
end


function constraint_on_off_des_pipe_head(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    # Get the design pipe from the specified index.
    des_pipe = ref(wm, nw, :des_pipe, a)

    # Initialize :on_off_des_pipe_head constraint dictionary.
    _initialize_con_dict(wm, :on_off_des_pipe_head, nw = nw, is_array = true)
    con(wm, nw, :on_off_des_pipe_head)[a] = Array{JuMP.ConstraintRef}([])

    # Apply the :on_off_des_pipe_head constraints.
    constraint_on_off_des_pipe_head(wm, nw, a, des_pipe["node_fr"], des_pipe["node_to"])
end


function constraint_on_off_des_pipe_head_loss(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    # Get the design pipe from the specified index.
    des_pipe = ref(wm, nw, :des_pipe, a)

    # Compute metadata associated with the design pipe.
    head_loss = wm.ref[:it][wm_it_sym][:head_loss]
    viscosity = wm.ref[:it][wm_it_sym][:viscosity]

    exponent = _get_exponent_from_head_loss_form(head_loss)

    wm_data = get_wm_data(wm.data)
    base_length = get(wm_data, "base_length", 1.0)
    base_mass = get(wm_data, "base_mass", 1.0)
    base_time = get(wm_data, "base_time", 1.0)

    r = _calc_pipe_resistance(
        des_pipe,
        head_loss,
        viscosity,
        base_length,
        base_mass,
        base_time,
    )
    q_max_reverse = min(get(des_pipe, "flow_max_reverse", 0.0), 0.0)
    q_min_forward = max(get(des_pipe, "flow_min_forward", 0.0), 0.0)

    # Initialize :on_off_des_pipe_head_loss constraint dictionary.
    _initialize_con_dict(wm, :on_off_des_pipe_head_loss, nw = nw, is_array = true)
    con(wm, nw, :on_off_des_pipe_head_loss)[a] = Array{JuMP.ConstraintRef}([])

    # Apply the :on_off_des_pipe_head_loss constraints.
    constraint_on_off_des_pipe_head_loss(
        wm,
        nw,
        a,
        des_pipe["node_fr"],
        des_pipe["node_to"],
        exponent,
        des_pipe["length"],
        r,
        q_max_reverse,
        q_min_forward,
    )
end


### Pump Constraints ###
function constraint_on_off_pump_flow(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr, node_to = ref(wm, nw, :pump, a)["node_fr"], ref(wm, nw, :pump, a)["node_to"]
    wm_data = get_wm_data(wm.data)
    flow_transform = _calc_flow_per_unit_transform(wm_data)
    q_min_forward = max(
        get(ref(wm, nw, :pump, a), "flow_min_forward", flow_transform(_FLOW_MIN)),
        flow_transform(_FLOW_MIN),
    )

    _initialize_con_dict(wm, :on_off_pump_flow, nw = nw, is_array = true)
    con(wm, nw, :on_off_pump_flow)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_pump_flow(wm, nw, a, q_min_forward)
end


function constraint_on_off_pump_head(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr, node_to = ref(wm, nw, :pump, a)["node_fr"], ref(wm, nw, :pump, a)["node_to"]

    _initialize_con_dict(wm, :on_off_pump_head, nw = nw, is_array = true)
    con(wm, nw, :on_off_pump_head)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_pump_head(wm, nw, a, node_fr, node_to)
end


function constraint_on_off_pump_head_gain(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr, node_to = ref(wm, nw, :pump, a)["node_fr"], ref(wm, nw, :pump, a)["node_to"]

    wm_data = get_wm_data(wm.data)
    flow_transform = _calc_flow_per_unit_transform(wm_data)
    q_min_forward = max(
        get(ref(wm, nw, :pump, a), "flow_min_forward", flow_transform(_FLOW_MIN)),
        flow_transform(_FLOW_MIN),
    )

    if ref(wm, nw, :pump, a, "pump_type") == PUMP_EPANET && isa(wm, AbstractNonlinearModel)
        message = "PUMP_EPANET head curves are not currently supported for nonlinear models."
        Memento.error(_LOGGER, message)
    end

    _initialize_con_dict(wm, :on_off_pump_head_gain, nw = nw, is_array = true)
    con(wm, nw, :on_off_pump_head_gain)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_pump_head_gain(wm, nw, a, node_fr, node_to, q_min_forward)
end


function constraint_on_off_pump_power(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    wm_data = get_wm_data(wm.data)
    flow_transform = _calc_flow_per_unit_transform(wm_data)
    q_min_forward = max(
        get(ref(wm, nw, :pump, a), "flow_min_forward", flow_transform(_FLOW_MIN)),
        flow_transform(_FLOW_MIN),
    )

    q_min_forward = min(q_min_forward, ref(wm, nw, :pump, a, "flow_max"))

    _initialize_con_dict(wm, :on_off_pump_power, nw = nw, is_array = true)
    con(wm, nw, :on_off_pump_power)[a] = Array{JuMP.ConstraintRef}([])

    if ref(wm, nw, :pump, a)["pump_type"] in [PUMP_QUADRATIC, PUMP_EPANET]
        # println("here")
        constraint_on_off_pump_power(wm, nw, a, q_min_forward)
    elseif ref(wm, nw, :pump, a)["pump_type"] == PUMP_BEST_EFFICIENCY_POINT
        wm_data = get_wm_data(wm.data)
        density = _calc_scaled_density(wm_data)
        gravity = _calc_scaled_gravity(wm_data)
        constraint_on_off_pump_power_best_efficiency(
            wm,
            nw,
            a,
            density,
            gravity,
            q_min_forward,
        )
    elseif ref(wm, nw, :pump, a)["pump_type"] == PUMP_LINEAR_POWER
        # Ensure that the required keys for modeling pump power exist.
        @assert haskey(ref(wm, nw, :pump, a), "power_fixed")
        @assert haskey(ref(wm, nw, :pump, a), "power_per_unit_flow")

        # Obtain the required data for modeling pump power linearly.
        power_fixed = ref(wm, nw, :pump, a)["power_fixed"]
        power_variable = ref(wm, nw, :pump, a)["power_per_unit_flow"]

        # Add the custom (linear) pump power constraint using the above.
        constraint_on_off_pump_power_custom(wm, nw, a, power_fixed, power_variable)
    end
end


function constraint_on_off_pump_group(
    wm::AbstractWaterModel,
    k::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    pump_indices = ref(wm, nw, :pump_group, k, "pump_indices")
    _initialize_con_dict(wm, :on_off_pump_group, nw = nw, is_array = true)
    con(wm, nw, :on_off_pump_group)[k] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_pump_group(wm, nw, k, Set{Int}(pump_indices))
end


function constraint_on_off_pump_switch(
    wm::AbstractWaterModel,
    a::Int,
    network_ids::Vector{Int};
    kwargs...,
)
    _initialize_con_dict(wm, :on_off_pump_switch, nw = network_ids[end], is_array = true)
    con(wm, network_ids[end], :on_off_pump_switch)[a] = Array{JuMP.ConstraintRef}([])
    max_switches = get(ref(wm, network_ids[end], :pump, a), "max_switches", 6)
    constraint_on_off_pump_switch(wm, a, network_ids, max_switches)
end


function constraint_pump_switch_on(
    wm::AbstractWaterModel,
    a::Int,
    n_1::Int,
    n_2::Int;
    kwargs...,
)
    _initialize_con_dict(wm, :pump_switch_on, nw = n_2, is_array = true)
    con(wm, n_2, :pump_switch_on)[a] = Array{JuMP.ConstraintRef}([])

    network_ids = sort(collect(nw_ids(wm)))
    min_active_time = get(ref(wm, n_2, :pump, a), "min_active_time", 3600.0)
    nw_end = n_2 + Int(floor(min_active_time / ref(wm, n_1, :time_step)))
    nws_active = Vector{Int}(collect(n_2:1:min(network_ids[end-1], nw_end)))
    constraint_pump_switch_on(wm, a, n_1, n_2, nws_active)
end


function constraint_pump_switch_off(
    wm::AbstractWaterModel,
    a::Int,
    n_1::Int,
    n_2::Int;
    kwargs...,
)
    _initialize_con_dict(wm, :pump_switch_off, nw = n_2, is_array = true)
    con(wm, n_2, :pump_switch_off)[a] = Array{JuMP.ConstraintRef}([])

    network_ids = sort(collect(nw_ids(wm)))
    min_inactive_time = get(ref(wm, n_2, :pump, a), "min_inactive_time", 1800.0)
    nw_end = n_2 + Int(floor(min_inactive_time / ref(wm, n_1, :time_step)))
    nws_inactive = Vector{Int}(collect(n_2:1:min(network_ids[end-1], nw_end)))
    constraint_pump_switch_off(wm, a, n_1, n_2, nws_inactive)
end

### Network expansion pump Constraints ###

function constraint_on_off_pump_flow_ne(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr, node_to = ref(wm, nw, :ne_pump, a)["node_fr"], ref(wm, nw, :ne_pump, a)["node_to"]
    wm_data = get_wm_data(wm.data)
    flow_transform = _calc_flow_per_unit_transform(wm_data)
    q_min_forward = max(
        get(ref(wm, nw, :ne_pump, a), "flow_min_forward", flow_transform(_FLOW_MIN)),
        flow_transform(_FLOW_MIN),
    )

    _initialize_con_dict(wm, :on_off_pump_flow_ne, nw = nw, is_array = true)
    con(wm, nw, :on_off_pump_flow_ne)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_pump_flow_ne(wm, nw, a, q_min_forward)
end

function constraint_on_off_pump_build_ne(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)

# Gather build and status variables
    x, z = var(wm, nw, :x_ne_pump, a), var(wm, nw, :z_ne_pump, a)
#
    _initialize_con_dict(wm, :on_off_pump_build_ne, nw = nw, is_array = true)
    con(wm, nw, :on_off_pump_build_ne)[a] = Array{JuMP.ConstraintRef}([])
    c = JuMP.@constraint(wm.model, z <= x)

    #Append the :on_off_pump_buld_ne constraint arrat
    append!(con(wm, nw, :on_off_pump_build_ne)[a], [c])
end

# function constraint_on_off_pump_build_ne(
#     wm::AbstractWaterModel,
#     a::Int,
#     nw::Int = nw_id_default,
#     kwargs...,
#     )
#
#     _initialize_con_dict(wm, :on_off_pump_buld_ne, nw = nw, is_array = true)
#     con(wm, nw, :on_off_pump_buld_ne)[a] = Array{JuMP.ConstraintRef}([])
#     constraint_on_off_pump_build_ne(wm, nw, a)
# end


# function constraint_on_off_pump_build_ne(
#     wm::AbstractWaterModel,
#     a::Int,
#     nw::Int = nw_id_default,
#     kwargs...,
#     )
#     println("testing $nw")
#     #Gather build and status variables
#     x, z = var(wm, nw, :x_ne_pump, a), var(wm, nw, :z_ne_pump, a)
#
#     con(wm, nw, :on_off_pump_build_ne)[a] = Array{JuMP.ConstraintRef}([])
#     c = JuMP.@constraint(wm.model, z <= x)
#
#     #Append the :on_off_pump_buld_ne constraint arrat
#     append!(con(wm, n, :on_off_pump_buld_ne)[a], [c])
# end

function constraint_on_off_pump_head_ne(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr, node_to = ref(wm, nw, :ne_pump, a)["node_fr"], ref(wm, nw, :ne_pump, a)["node_to"]

    _initialize_con_dict(wm, :on_off_pump_head_ne, nw = nw, is_array = true)
    con(wm, nw, :on_off_pump_head_ne)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_pump_head_ne(wm, nw, a, node_fr, node_to)
end

function constraint_on_off_pump_head_gain_ne(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr, node_to = ref(wm, nw, :ne_pump, a)["node_fr"], ref(wm, nw, :ne_pump, a)["node_to"]

    wm_data = get_wm_data(wm.data)
    flow_transform = _calc_flow_per_unit_transform(wm_data)
    q_min_forward = max(
        get(ref(wm, nw, :ne_pump, a), "flow_min_forward", flow_transform(_FLOW_MIN)),
        flow_transform(_FLOW_MIN),
    )

    if ref(wm, nw, :ne_pump, a, "pump_type") == PUMP_EPANET && isa(wm, AbstractNonlinearModel)
        message = "NE PUMP_EPANET head curves are not currently supported for nonlinear models."
        Memento.error(_LOGGER, message)
    end

    _initialize_con_dict(wm, :on_off_pump_head_gain_ne, nw = nw, is_array = true)
    con(wm, nw, :on_off_pump_head_gain_ne)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_pump_head_gain_ne(wm, nw, a, node_fr, node_to, q_min_forward)
end

function constraint_on_off_pump_power_ne(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    wm_data = get_wm_data(wm.data)
    flow_transform = _calc_flow_per_unit_transform(wm_data)
    q_min_forward = max(
        get(ref(wm, nw, :ne_pump, a), "flow_min_forward", flow_transform(_FLOW_MIN)),
        flow_transform(_FLOW_MIN),
    )

    q_min_forward = min(q_min_forward, ref(wm, nw, :ne_pump, a, "flow_max"))

    _initialize_con_dict(wm, :on_off_pump_power_ne, nw = nw, is_array = true)
    con(wm, nw, :on_off_pump_power_ne)[a] = Array{JuMP.ConstraintRef}([])

    if ref(wm, nw, :ne_pump, a)["pump_type"] in [PUMP_QUADRATIC, PUMP_EPANET]
        constraint_on_off_pump_power_ne(wm, nw, a, q_min_forward)
    elseif ref(wm, nw, :ne_pump, a)["pump_type"] == PUMP_BEST_EFFICIENCY_POINT
        wm_data = get_wm_data(wm.data)
        density = _calc_scaled_density(wm_data)
        gravity = _calc_scaled_gravity(wm_data)
        constraint_on_off_pump_power_best_efficiency_ne(
            wm,
            nw,
            a,
            density,
            gravity,
            q_min_forward,
        )
    elseif ref(wm, nw, :ne_pump, a)["pump_type"] == PUMP_LINEAR_POWER
        # Ensure that the required keys for modeling expansion pump power exist.
        @assert haskey(ref(wm, nw, :ne_pump, a), "power_fixed")
        @assert haskey(ref(wm, nw, :ne_pump, a), "power_per_unit_flow")

        # Obtain the required data for modeling expansion pump power linearly.
        power_fixed = ref(wm, nw, :ne_pump, a)["power_fixed"]
        power_variable = ref(wm, nw, :ne_pump, a)["power_per_unit_flow"]

        # Add the custom (linear) pump power constraint using the above.
        constraint_on_off_pump_power_custom_ne(wm, nw, a, power_fixed, power_variable)
    end
end

function constraint_on_off_pump_group_ne(
    wm::AbstractWaterModel,
    k::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    pump_indices = ref(wm, nw, :ne_pump_group, k, "pump_indices")
    _initialize_con_dict(wm, :on_off_pump_group_ne, nw = nw, is_array = true)
    con(wm, nw, :on_off_pump_group_ne)[k] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_pump_group_ne(wm, nw, k, Set{Int}(pump_indices))
end

function constraint_on_off_pump_switch_ne(
    wm::AbstractWaterModel,
    a::Int,
    network_ids::Vector{Int};
    kwargs...,
)
    _initialize_con_dict(wm, :on_off_pump_switch_ne, nw = network_ids[end], is_array = true)
    con(wm, network_ids[end], :on_off_pump_switch_ne)[a] = Array{JuMP.ConstraintRef}([])
    max_switches = get(ref(wm, network_ids[end], :ne_pump, a), "max_switches", 6)
    constraint_on_off_pump_switch_ne(wm, a, network_ids, max_switches)
end

function constraint_pump_switch_on_ne(
    wm::AbstractWaterModel,
    a::Int,
    n_1::Int,
    n_2::Int;
    kwargs...,
)
    _initialize_con_dict(wm, :ne_pump_switch_on, nw = n_2, is_array = true)
    con(wm, n_2, :ne_pump_switch_on)[a] = Array{JuMP.ConstraintRef}([])

    network_ids = sort(collect(nw_ids(wm)))
    min_active_time = get(ref(wm, n_2, :ne_pump, a), "min_active_time", 3600.0)
    nw_end = n_2 + Int(floor(min_active_time / ref(wm, n_1, :time_step)))
    nws_active = Vector{Int}(collect(n_2:1:min(network_ids[end-1], nw_end)))
    constraint_pump_switch_on_ne(wm, a, n_1, n_2, nws_active)
end

function constraint_pump_switch_off_ne(
    wm::AbstractWaterModel,
    a::Int,
    n_1::Int,
    n_2::Int;
    kwargs...,
)
    _initialize_con_dict(wm, :pump_switch_off_ne, nw = n_2, is_array = true)
    con(wm, n_2, :pump_switch_off_ne)[a] = Array{JuMP.ConstraintRef}([])

    network_ids = sort(collect(nw_ids(wm)))
    min_inactive_time = get(ref(wm, n_2, :ne_pump, a), "min_inactive_time", 1800.0)
    nw_end = n_2 + Int(floor(min_inactive_time / ref(wm, n_1, :time_step)))
    nws_inactive = Vector{Int}(collect(n_2:1:min(network_ids[end-1], nw_end)))
    constraint_pump_switch_off_ne(wm, a, n_1, n_2, nws_inactive)
end


### Short Pipe Constraints ###
function constraint_short_pipe_flow(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    short_pipe = ref(wm, nw, :short_pipe, a)
    q_max_reverse = min(get(short_pipe, "flow_max_reverse", 0.0), short_pipe["flow_max"])
    q_min_forward = max(get(short_pipe, "flow_min_forward", 0.0), short_pipe["flow_min"])

    _initialize_con_dict(wm, :short_pipe_flow, nw = nw, is_array = true)
    con(wm, nw, :short_pipe_flow)[a] = Array{JuMP.ConstraintRef}([])
    constraint_short_pipe_flow(wm, nw, a, q_max_reverse, q_min_forward)
end


function constraint_short_pipe_head(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr = ref(wm, nw, :short_pipe, a)["node_fr"]
    node_to = ref(wm, nw, :short_pipe, a)["node_to"]

    _initialize_con_dict(wm, :short_pipe_head, nw = nw, is_array = true)
    con(wm, nw, :short_pipe_head)[a] = Array{JuMP.ConstraintRef}([])
    constraint_short_pipe_head(wm, nw, a, node_fr, node_to)
end


### Network Expansion Short Pipe Constraints ###
function constraint_short_pipe_flow_ne(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    ne_short_pipe = ref(wm, nw, :ne_short_pipe, a)
    q_max_reverse = min(get(ne_short_pipe, "flow_max_reverse", 0.0), ne_short_pipe["flow_max"])
    q_min_forward = max(get(ne_short_pipe, "flow_min_forward", 0.0), ne_short_pipe["flow_min"])

    _initialize_con_dict(wm, :short_pipe_flow_ne, nw = nw, is_array = true)
    con(wm, nw, :short_pipe_flow_ne)[a] = Array{JuMP.ConstraintRef}([])
    constraint_short_pipe_flow_ne(wm, nw, a, q_max_reverse, q_min_forward)
end


function constraint_short_pipe_head_ne(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr = ref(wm, nw, :ne_short_pipe, a)["node_fr"]
    node_to = ref(wm, nw, :ne_short_pipe, a)["node_to"]

    _initialize_con_dict(wm, :short_pipe_head_ne, nw = nw, is_array = true)
    con(wm, nw, :short_pipe_head_ne)[a] = Array{JuMP.ConstraintRef}([])
    constraint_short_pipe_head_ne(wm, nw, a, node_fr, node_to)
end


### Valve Constraints ###
function constraint_on_off_valve_flow(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    valve = ref(wm, nw, :valve, a)
    node_fr, node_to = valve["node_fr"], valve["node_to"]
    q_max_reverse = min(get(valve, "flow_max_reverse", 0.0), valve["flow_max"])
    q_min_forward = max(get(valve, "flow_min_forward", 0.0), valve["flow_min"])

    _initialize_con_dict(wm, :on_off_valve_flow, nw = nw, is_array = true)
    con(wm, nw, :on_off_valve_flow)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_valve_flow(wm, nw, a, q_max_reverse, q_min_forward)
end


function constraint_on_off_valve_head(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr, node_to = ref(wm, nw, :valve, a)["node_fr"], ref(wm, nw, :valve, a)["node_to"]

    _initialize_con_dict(wm, :on_off_valve_head, nw = nw, is_array = true)
    con(wm, nw, :on_off_valve_head)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_valve_head(wm, nw, a, node_fr, node_to)
end


### Regulator Constraints ###
function constraint_on_off_regulator_flow(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    wm_data = get_wm_data(wm.data)
    flow_transform = _calc_flow_per_unit_transform(wm_data)
    q_min_forward = max(
        get(ref(wm, nw, :regulator, a), "flow_min_forward", flow_transform(_FLOW_MIN)),
        flow_transform(_FLOW_MIN),
    )

    _initialize_con_dict(wm, :on_off_regulator_flow, nw = nw, is_array = true)
    con(wm, nw, :on_off_regulator_flow)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_regulator_flow(wm, nw, a, q_min_forward)
end


function constraint_on_off_regulator_head(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr = ref(wm, nw, :regulator, a)["node_fr"]
    node_to = ref(wm, nw, :regulator, a)["node_to"]
    elevation = ref(wm, nw, :node, node_to)["elevation"]
    head_setting = elevation + ref(wm, nw, :regulator, a, "setting")

    _initialize_con_dict(wm, :on_off_regulator_head, nw = nw, is_array = true)
    con(wm, nw, :on_off_regulator_head)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_regulator_head(wm, nw, a, node_fr, node_to, head_setting)
end
