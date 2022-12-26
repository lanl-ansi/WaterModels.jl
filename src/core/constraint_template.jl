# Constraint templates help simplify data wrangling across multiple water network
# optimization formulations by providing an abstraction layer between the network data and
# network constraint definitions. Each constraint template's job is to extract the required
# parameters from a given network data structure and pass the data as named arguments to
# the water network optimization constraints. Here, constraint templates should always be
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

"""
    constraint_flow_conservation(
        wm::AbstractWaterModel,
        i::Int;
        nw::Int=nw_id_default
    )

Constraint template to add constraints that ensure volumetric flow rate (and
thus mass) is conserved at node `i` and subnetwork (or time) index `nw` in the
network. Here, `wm` is the WaterModels object.
"""
function constraint_flow_conservation(
    wm::AbstractWaterModel,
    i::Int;
    nw::Int = nw_id_default,
)
    # Collect various indices for edge-type components connected to node `i`.
    pipe_fr, pipe_to = _get_from_and_to_components(wm, i, :pipe; nw = nw)
    des_pipe_fr, des_pipe_to = _get_from_and_to_components(wm, i, :des_pipe; nw = nw)
    pump_fr, pump_to = _get_from_and_to_components(wm, i, :pump; nw = nw)
    regulator_fr, regulator_to = _get_from_and_to_components(wm, i, :regulator; nw = nw)
    short_pipe_fr, short_pipe_to = _get_from_and_to_components(wm, i, :short_pipe; nw = nw)
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
        regulator_fr,
        regulator_to,
        short_pipe_fr,
        short_pipe_to,
        valve_fr,
        valve_to,
        reservoirs,
        tanks,
        dispatchable_demands,
        net_fixed_demand,
    )
end


"""
    constraint_node_directionality(
        wm::AbstractWaterModel,
        i::Int;
        nw::Int=nw_id_default
    )

Constraint template to add direction-based constraints
([`constraint_sink_directionality`](@ref),
[`constraint_source_directionality`](@ref), or
[`constraint_intermediate_directionality`](@ref)) when appropriate. Here, `wm`
is the WaterModels object, `i` is the index of the node for which the
constraints will be added, if applicable, and `nw` is the subnetwork (or time)
index.
"""
function constraint_node_directionality(
    wm::AbstractWaterModel,
    i::Int;
    nw::Int = nw_id_default,
)
    # Collect various indices for edge-type components connected to node `i`.
    pipe_fr, pipe_to = _get_from_and_to_components(wm, i, :pipe; nw = nw)
    des_pipe_fr, des_pipe_to = _get_from_and_to_components(wm, i, :des_pipe; nw = nw)
    pump_fr, pump_to = _get_from_and_to_components(wm, i, :pump; nw = nw)
    regulator_fr, regulator_to = _get_from_and_to_components(wm, i, :regulator; nw = nw)
    short_pipe_fr, short_pipe_to = _get_from_and_to_components(wm, i, :short_pipe; nw = nw)
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
        length(regulator_to) +
        length(short_pipe_to) +
        length(valve_to)

    # Get the out degree of node `i`.
    out_length =
        length(pipe_fr) +
        length(des_pipe_fr) +
        length(pump_fr) +
        length(regulator_fr) +
        length(short_pipe_fr) +
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
            regulator_fr,
            regulator_to,
            short_pipe_fr,
            short_pipe_to,
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
            regulator_fr,
            regulator_to,
            short_pipe_fr,
            short_pipe_to,
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
            regulator_fr,
            regulator_to,
            short_pipe_fr,
            short_pipe_to,
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
            regulator_fr,
            regulator_to,
            short_pipe_fr,
            short_pipe_to,
            valve_fr,
            valve_to,
        )
    end
end


### Tank Constraints ###

"""
    constraint_tank_volume(
        wm::AbstractWaterModel,
        i::Int;
        nw::Int=nw_id_default
    )

Constraint template to add [`constraint_tank_volume_fixed`](@ref) constraints
when a tank is not dispatchable, usually at the first time index of a problem
to fix the tank's initial volume. Here, `wm` is the WaterModels object, `i` is
the index of the tank for which the constraints will be added, if applicable,
and `nw` is the subnetwork (or time) index at which the volume will be fixed.
"""
function constraint_tank_volume(wm::AbstractWaterModel, i::Int; nw::Int = nw_id_default)
    # Only set the tank state if the tank is nondispatchable.
    if !ref(wm, nw, :tank, i)["dispatchable"]
        tank = ref(wm, nw, :tank, i)
        initial_level = tank["init_level"]
        surface_area = 0.25 * pi * tank["diameter"]^2

        V_initial = surface_area * initial_level
        time_step = ref(wm, nw, :time_step)
        V_min = max(tank["min_vol"], surface_area * tank["min_level"])
        V_max = surface_area * tank["max_level"]

        # Update the nodal elevation data.
        node = ref(wm, nw, :node, tank["node"])
        head = node["elevation"] + initial_level
        node["head_nominal"] = node["head_min"] = node["head_max"] = head

        # Apply the tank volume constraint at the specified time step.
        _initialize_con_dict(wm, :tank_volume; nw = nw, is_array = true)
        con(wm, nw, :tank_volume)[i] = Array{JuMP.ConstraintRef,1}([])
        constraint_tank_volume_fixed(wm, nw, i, V_initial, time_step, V_min, V_max)
    end
end


"""
    constraint_tank_volume(
        wm::AbstractWaterModel,
        i::Int,
        nw_1::Int,
        nw_2::Int
    )

Constraint template to add [`constraint_tank_volume`](@ref) constraints that
integrate the volume of a tank forward in time. Here, `wm` is the WaterModels
object, `i` is the index of the tank for which the constraints will be added,
`nw_1` is the first time index considered in the temporal integration, and
`nw_2` is the adjacent (second) time index considered in the integration.
"""
function constraint_tank_volume(wm::AbstractWaterModel, i::Int, nw_1::Int, nw_2::Int)
    # Only apply the constraint if the tank exists in both subnetworks.
    if haskey(ref(wm, nw_1, :tank), i) && haskey(ref(wm, nw_2, :tank), i)
        # Get the tank reference within each of the subnetworks.
        tank_nw_1, tank_nw_2 = ref(wm, nw_1, :tank, i), ref(wm, nw_2, :tank, i)

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

"""
    constraint_pipe_flow(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_pipe_flow`](@ref) constraints that
limit the volumetric flow rate across a pipe. Here, `wm` is the WaterModels
object, `a` is the index of the pipe for which flow will be limited, and `nw`
is the subnetwork (or time) index that is considered.
"""
function constraint_pipe_flow(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    q_max_reverse = min(get(ref(wm, nw, :pipe, a), "flow_max_reverse", 0.0), 0.0)
    q_min_forward = max(get(ref(wm, nw, :pipe, a), "flow_min_forward", 0.0), 0.0)

    _initialize_con_dict(wm, :pipe_flow, nw = nw, is_array = true)
    con(wm, nw, :pipe_flow)[a] = Array{JuMP.ConstraintRef}([])
    constraint_pipe_flow(wm, nw, a, q_max_reverse, q_min_forward)
end


"""
    constraint_pipe_head(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_pipe_head`](@ref) constraints that
limit and establish relationships among head difference and head variables.
Here, `wm` is the WaterModels object, `a` is the index of the pipe, and `nw` is
the subnetwork (or time) index that is considered.
"""
function constraint_pipe_head(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr, node_to = ref(wm, nw, :pipe, a)["node_fr"], ref(wm, nw, :pipe, a)["node_to"]
    _initialize_con_dict(wm, :pipe_head, nw = nw, is_array = true)
    con(wm, nw, :pipe_head)[a] = Array{JuMP.ConstraintRef}([])
    constraint_pipe_head(wm, nw, a, node_fr, node_to)
end


"""
    constraint_pipe_head_loss(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_pipe_head_loss`](@ref) constraints
that model head loss relationships along a pipe. Here, `wm` is the WaterModels
object, `a` is the index of the pipe, and `nw` is the subnetwork (or time)
index that is considered.
"""
function constraint_pipe_head_loss(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr, node_to = ref(wm, nw, :pipe, a)["node_fr"], ref(wm, nw, :pipe, a)["node_to"]
    exponent, L = ref(wm, nw, :alpha), ref(wm, nw, :pipe, a)["length"]

    wm_data = get_wm_data(wm.data)
    head_loss, viscosity = wm_data["head_loss"], wm_data["viscosity"]
    base_length = wm_data["per_unit"] ? wm_data["base_length"] : 1.0
    base_mass = wm_data["per_unit"] ? wm_data["base_mass"] : 1.0
    base_time = wm_data["per_unit"] ? wm_data["base_time"] : 1.0

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

"""
    constraint_des_pipe_flow(
        wm::AbstractWaterModel,
        k::Int,
        node_fr::Int,
        node_to::Int;
        nw::Int=nw_id_default,
        kwargs...,
    )

Constraint template to add constraints that ensure flow variables for design
pipes along an arc are equally directed (for direction-based formulations).
Here, `wm` is the WaterModels object, `k` is the index of the arc that connects
the two common nodes of each design pipe, `node_fr` is the index of the tail
node of the arc that models each design pipe, `node_to` is the index of the
head node of the arc that models each design pipe, `nw` is the index of a
subnetwork within a multinetwork, and `des_pipes` is the vector of design pipe
indices for design pipes that reside along the same arc `k`.
"""
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


"""
    constraint_des_pipe_head(
        wm::AbstractWaterModel,
        k::Int,
        node_fr::Int,
        node_to::Int;
        nw::Int=nw_id_default,
        kwargs...,
    )

Constraint template to add constraints that ensure head difference variables
for design pipes along an arc sum to the actual head difference along that
arc. Here, `wm` is the WaterModels object, `k` is the index of the arc that
connects the two common nodes of each design pipe, `node_fr` is the index of
the tail node of the arc that models each design pipe, `node_to` is the index
of the head node of the arc that models each design pipe, `nw` is the index of
a subnetwork within a multinetwork, and `des_pipes` is the vector of design
pipe indices for design pipes that reside along the same arc `k`.
"""
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


"""
    constraint_des_pipe_selection(
        wm::AbstractWaterModel,
        k::Int,
        node_fr::Int,
        node_to::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_des_pipe_selection`](@ref) constraints
that enforce the selection of only one design pipe to be constructed along a
given arc. Here, `wm` is the WaterModels object, `k` is the index of the arc
that connects the two common nodes of each design pipe, `node_fr` is the index
of the tail node of the arc that models each design pipe, `node_to` is the
index of the head node of the arc that models each design pipe, and `nw` is the
index of a subnetwork within a multinetwork.
"""
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


"""
    constraint_on_off_des_pipe_flow(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_on_off_des_pipe_flow`](@ref)
constraints that limit the amount of flow traveling across a design pipe. Here,
`wm` is the WaterModels object, `a` is the index of the design pipe, and `nw`
is the index of a subnetwork within a multinetwork.
"""
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


"""
    constraint_on_off_des_pipe_head(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_on_off_des_pipe_head`](@ref)
constraints that limit the head differences across a design pipe. Here, `wm` is
the WaterModels object, `a` is the index of the design pipe, and `nw` is the
index of a subnetwork within a multinetwork.
"""
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


"""
    constraint_on_off_des_pipe_head_loss(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_on_off_des_pipe_head_loss`](@ref)
constraints that model the head losses across a design pipe. Here, `wm` is the
WaterModels object, `a` is the index of the design pipe, and `nw` is the index
of a subnetwork within a multinetwork.
"""
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
    base_length = wm_data["per_unit"] ? wm_data["base_length"] : 1.0
    base_mass = wm_data["per_unit"] ? wm_data["base_mass"] : 1.0
    base_time = wm_data["per_unit"] ? wm_data["base_time"] : 1.0

    r = _calc_pipe_resistance(des_pipe, head_loss, viscosity, base_length, base_mass, base_time)
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

"""
    constraint_on_off_pump_flow(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_on_off_pump_flow`](@ref) constraints,
which restrict the amount of flow transported through a pump based on its
operating status. Here, `wm` is the WaterModels object, `a` is the index of the
pump, and `nw` is the index of a subnetwork within a multinetwork.
"""
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


"""
    constraint_on_off_pump_flow(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_on_off_pump_head`](@ref) constraints,
which disjunctively limit the head difference between nodes connected by the
pump and, if operating, ensures the head difference between nodes is equal to
the head gain, constrained by ([`constraint_on_off_pump_head_gain`](@ref).
Here, `wm` is the WaterModels object, `a` is the index of the pump, and `nw` is
the index of a subnetwork within a multinetwork.
"""
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


"""
    constraint_on_off_pump_head_gain(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_on_off_pump_head_gain`](@ref)
constraints, which, if operating, limit the pump's head gain as a function of
flow rate. Here, `wm` is the WaterModels object, `a` is the index of the pump,
and `nw` is the index of a subnetwork within a multinetwork.
"""
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

    if ref(wm, nw, :pump, a)["head_curve_form"] == PUMP_EPANET &&
       isa(wm, AbstractNonlinearModel)
        message = "PUMP_EPANET head curves are not currently supported for nonlinear models."
        Memento.error(_LOGGER, message)
    end

    _initialize_con_dict(wm, :on_off_pump_head_gain, nw = nw, is_array = true)
    con(wm, nw, :on_off_pump_head_gain)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_pump_head_gain(wm, nw, a, node_fr, node_to, q_min_forward)
end


"""
    constraint_on_off_pump_power(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_on_off_pump_power`](@ref) constraints,
which, if operating, model the pump's power according to certain assumptions.
Here, `wm` is the WaterModels object, `a` is the index of the pump, and `nw` is
the index of a subnetwork within a multinetwork.
"""
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

    _initialize_con_dict(wm, :on_off_pump_power, nw = nw, is_array = true)
    con(wm, nw, :on_off_pump_power)[a] = Array{JuMP.ConstraintRef}([])

    if ref(wm, nw, :pump, a)["head_curve_form"] in [PUMP_QUADRATIC, PUMP_EPANET]
        constraint_on_off_pump_power(wm, nw, a, q_min_forward)
    elseif ref(wm, nw, :pump, a)["head_curve_form"] == PUMP_BEST_EFFICIENCY_POINT
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
    elseif ref(wm, nw, :pump, a)["head_curve_form"] == PUMP_LINEAR_POWER
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


"""
    constraint_on_off_pump_power(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_on_off_pump_group`](@ref) constraints,
which impose symmetry-breaking lexicographic sorting of pump activation
statuses on groups of identical pumps operating in parallel along the same arc
of the network. Here, `wm` is the WaterModels object, `k` is the index of the
pump group, and `nw` is the index of a subnetwork within a multinetwork.
"""
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


"""
    constraint_on_off_pump_switch(
        wm::AbstractWaterModel,
        a::Int;
        network_ids::Array{Int,1};
        kwargs...
    )

Constraint template to add [`constraint_on_off_pump_switch`](@ref) constraint,
which limits the number of times a pump can be switched from off to on in a
multiperiod pump scheduling problem. Here, `wm` is the WaterModels object, `a`
is the index of the pump, and `network_ids` are the network (time) indices used
in the summation that limits the number of switches.
"""
function constraint_on_off_pump_switch(
    wm::AbstractWaterModel,
    a::Int,
    network_ids::Array{Int,1};
    kwargs...,
)
    _initialize_con_dict(wm, :on_off_pump_switch, nw = network_ids[end], is_array = true)
    con(wm, network_ids[end], :on_off_pump_switch)[a] = Array{JuMP.ConstraintRef}([])
    max_switches = get(ref(wm, network_ids[end], :pump, a), "max_switches", 6)
    constraint_on_off_pump_switch(wm, a, network_ids, max_switches)
end


"""
    constraint_pump_switch_on(
        wm::AbstractWaterModel,
        a::Int,
        n_1::Int,
        n_2::Int;
        kwargs...
    )

Constraint template to add [`constraint_pump_switch_on`](@ref) constraints,
which model the switching of a pump from _off_ to _on_ and constrains its
operation, if switched on, to remain on for at least some minimum amount of
time. Here, `wm` is the WaterModels object; `a` is the index of the pump; `n_1`
is the first time index modeled by the constraint; and `n_2` is the adjacent
(next) time index modeled by the constraint, which permits limiting the change
in pump status (i.e., from off to on).
"""
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


"""
    function constraint_pump_switch_off(
        wm::AbstractWaterModel,
        a::Int,
        n_1::Int,
        n_2::Int;
        kwargs...
    )

Constraint template to add [`constraint_pump_switch_off`](@ref) constraints,
which model the switching of a pump from _on_ to _off_ and constrains non-
operation, if switched off, to remain off for at least some minimum amount of
time. Here, `wm` is the WaterModels object; `a` is the index of the pump; `n_1`
is the first time index modeled by the constraint; `n_2` is the adjacent (next)
time index modeled by the constraint, which permits limiting the change in pump
status (i.e., from on to off).
"""
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


### Short Pipe Constraints ###

"""
    constraint_short_pipe_flow(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_short_pipe_flow`](@ref) constraints
that limit the volumetric flow rate across a short pipe. Here, `wm` is the
WaterModels object, `a` is the index of the short pipe for which flow will be
limited, and `nw` is the subnetwork (or time) index that is considered.
"""
function constraint_short_pipe_flow(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    q_max_reverse = min(get(ref(wm, nw, :short_pipe, a), "flow_max_reverse", 0.0), 0.0)
    q_min_forward = max(get(ref(wm, nw, :short_pipe, a), "flow_min_forward", 0.0), 0.0)

    _initialize_con_dict(wm, :short_pipe_flow, nw = nw, is_array = true)
    con(wm, nw, :short_pipe_flow)[a] = Array{JuMP.ConstraintRef}([])
    constraint_short_pipe_flow(wm, nw, a, q_max_reverse, q_min_forward)
end


"""
    constraint_short_pipe_head(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_short_pipe_head`](@ref) constraints
that limit and establish relationships among head difference and head
variables. Here, `wm` is the WaterModels object, `a` is the index of the short
pipe, and `nw` is the subnetwork (or time) index that is considered.
"""
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


### Valve Constraints ###

"""
    constraint_on_off_valve_flow(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_on_off_valve_flow`](@ref) constraints
that limit the volumetric flow rate across a valve based on its operating
status. Here, `wm` is the WaterModels object, `a` is the index of the valve for
which flow will be limited, and `nw` is the subnetwork (or time) index that is
considered.
"""
function constraint_on_off_valve_flow(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr, node_to = ref(wm, nw, :valve, a)["node_fr"], ref(wm, nw, :valve, a)["node_to"]
    q_max_reverse = min(get(ref(wm, nw, :valve, a), "flow_max_reverse", 0.0), 0.0)
    q_min_forward = max(get(ref(wm, nw, :valve, a), "flow_min_forward", 0.0), 0.0)

    _initialize_con_dict(wm, :on_off_valve_flow, nw = nw, is_array = true)
    con(wm, nw, :on_off_valve_flow)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_valve_flow(wm, nw, a, q_max_reverse, q_min_forward)
end


"""
    constraint_on_off_valve_head(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_on_off_valve_head`](@ref) constraints
that limit and establish relationships among head variables based on the
operating status of the valve. Here, `wm` is the WaterModels object, `a` is the
index of the valve, and `nw` is the subnetwork (or time) index that is
considered.
"""
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

"""
    constraint_on_off_regulator_flow(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_on_off_regulator_flow`](@ref)
constraints that limit the volumetric flow rate across a regulator based on its
operating status. Here, `wm` is the WaterModels object, `a` is the index of the
valve for which flow will be limited, and `nw` is the subnetwork (or time)
index that is considered.
"""
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


"""
    constraint_on_off_regulator_head(
        wm::AbstractWaterModel,
        a::Int;
        nw::Int=nw_id_default,
        kwargs...
    )

Constraint template to add [`constraint_on_off_regulator_head`](@ref)
constraints that limit and establish relationships among head variables based
on the operating status of the regulator. Here, `wm` is the WaterModels object,
`a` is the index of the regulator, and `nw` is the subnetwork (or time) index
that is considered.
"""
function constraint_on_off_regulator_head(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr = ref(wm, nw, :regulator, a)["node_fr"]
    node_to = ref(wm, nw, :regulator, a)["node_to"]
    elevation = ref(wm, nw, :node, node_to)["elevation"]
    head_setting = elevation + ref(wm, nw, :regulator, a)["setting"]

    _initialize_con_dict(wm, :on_off_regulator_head, nw = nw, is_array = true)
    con(wm, nw, :on_off_regulator_head)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_regulator_head(wm, nw, a, node_fr, node_to, head_setting)
end
