# Constraint templates help simplify data wrangling across multiple Water Flow
# formulations by providing an abstraction layer between the network data and
# network constraint definitions. The constraint template's job is to extract
# the required parameters from a given network data structure and pass the data
# as named arguments to the Water Flow formulations.
#
# Constraint templates should always be defined over "AbstractWaterModel" and
# should never refer to model variables.

function _initialize_con_dict(wm::AbstractWaterModel, key::Symbol; nw::Int=wm.cnw, is_array::Bool=false)
    if !haskey(con(wm, nw), key)
        con(wm, nw)[key] = is_array ? Dict{Int, Array{JuMP.ConstraintRef}}() :
            Dict{Int, JuMP.ConstraintRef}()
    end
end


function _collect_comps_fr(wm::AbstractWaterModel, i::Int, sym::Symbol; nw::Int=wm.cnw)
    return collect(keys(filter(x -> x.second["node_fr"] == i, ref(wm, nw, sym))))
end


function _collect_comps_to(wm::AbstractWaterModel, i::Int, sym::Symbol; nw::Int=wm.cnw)
    return collect(keys(filter(x -> x.second["node_to"] == i, ref(wm, nw, sym))))
end


### Nodal Constraints ###
function constraint_flow_conservation(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    # Collect various indices for edge-type components connected to node `i`.
    pipe_fr = _collect_comps_fr(wm, i, :pipe; nw=nw)
    pipe_to = _collect_comps_to(wm, i, :pipe; nw=nw)
    des_pipe_fr = _collect_comps_fr(wm, i, :des_pipe; nw=nw)
    des_pipe_to = _collect_comps_to(wm, i, :des_pipe; nw=nw)
    pump_fr = _collect_comps_fr(wm, i, :pump; nw=nw)
    pump_to = _collect_comps_to(wm, i, :pump; nw=nw)
    regulator_fr = _collect_comps_fr(wm, i, :regulator; nw=nw)
    regulator_to = _collect_comps_to(wm, i, :regulator; nw=nw)
    short_pipe_fr = _collect_comps_fr(wm, i, :short_pipe; nw=nw)
    short_pipe_to = _collect_comps_to(wm, i, :short_pipe; nw=nw)
    valve_fr = _collect_comps_fr(wm, i, :valve; nw=nw)
    valve_to = _collect_comps_to(wm, i, :valve; nw=nw)

    # Collect various indices for node-type components connected to node `i`.
    reservoirs = ref(wm, nw, :node_reservoir, i) # Reservoirs attached to node `i`.
    tanks = ref(wm, nw, :node_tank, i) # Tanks attached to node `i`.
    demands = ref(wm, nw, :node_demand, i) # Demands attached to node `i`.

    # Sum the constant demands required at node `i`.
    nondispatchable_demands = filter(j -> j in ids(wm, nw, :nondispatchable_demand), demands)
    fixed_demands = [ref(wm, nw, :nondispatchable_demand, j)["flow_rate"] for j in nondispatchable_demands]
    net_fixed_demand = length(fixed_demands) > 0 ? sum(fixed_demands) : 0.0

    # Get the indices of dispatchable demands connected to node `i`.
    dispatchable_demands = filter(j -> j in ids(wm, nw, :dispatchable_demand), demands)

    # Initialize the flow conservation constraint dictionary entry.
    _initialize_con_dict(wm, :flow_conservation, nw=nw)

    # Add the flow conservation constraint.
    constraint_flow_conservation(
        wm, nw, i, pipe_fr, pipe_to, des_pipe_fr, des_pipe_to, pump_fr, pump_to,
        regulator_fr, regulator_to, short_pipe_fr, short_pipe_to, valve_fr, valve_to,
        reservoirs, tanks, dispatchable_demands, net_fixed_demand)
end


function constraint_node_directionality(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    # Collect various indices for edge-type components connected to node `i`.
    pipe_fr = _collect_comps_fr(wm, i, :pipe; nw=nw)
    pipe_to = _collect_comps_to(wm, i, :pipe; nw=nw)
    des_pipe_fr = _collect_comps_fr(wm, i, :des_pipe; nw=nw)
    des_pipe_to = _collect_comps_to(wm, i, :des_pipe; nw=nw)
    pump_fr = _collect_comps_fr(wm, i, :pump; nw=nw)
    pump_to = _collect_comps_to(wm, i, :pump; nw=nw)
    regulator_fr = _collect_comps_fr(wm, i, :regulator; nw=nw)
    regulator_to = _collect_comps_to(wm, i, :regulator; nw=nw)
    short_pipe_fr = _collect_comps_fr(wm, i, :short_pipe; nw=nw)
    short_pipe_to = _collect_comps_to(wm, i, :short_pipe; nw=nw)
    valve_fr = _collect_comps_fr(wm, i, :valve; nw=nw)
    valve_to = _collect_comps_to(wm, i, :valve; nw=nw)

    # Collect various indices for node-type components connected to node `i`.
    reservoirs = ref(wm, nw, :node_reservoir, i) # Reservoirs attached to node `i`.
    tanks = ref(wm, nw, :node_tank, i) # Tanks attached to node `i`.
    demands = ref(wm, nw, :node_demand, i) # Demands attached to node `i`.

    # Sum the constant demands required at node `i`.
    nondispatchable_demands = filter(j -> j in ids(wm, nw, :nondispatchable_demand), demands)
    fixed_demands = [ref(wm, nw, :nondispatchable_demand, j)["flow_rate"] for j in nondispatchable_demands]
    net_fixed_demand = length(fixed_demands) > 0 ? sum(fixed_demands) : 0.0

    # Get the number of nodal components attached to node `i`.
    num_components = length(demands) + length(tanks) + length(reservoirs)

    # Get the in degree of node `i`.
    in_length = length(pipe_to) + length(des_pipe_to) + length(pump_to) +
        length(regulator_to) + length(short_pipe_to) + length(valve_to)

    # Get the out degree of node `i`.
    out_length = length(pipe_fr) + length(des_pipe_fr) + length(pump_fr) +
        length(regulator_fr) + length(short_pipe_fr) + length(valve_fr)

    # Initialize the directionality constraint dictionary entry.
    _initialize_con_dict(wm, :node_directionality, nw=nw)

    # Check if node directionality constraints should be added.
    if num_components == 0 && in_length + out_length == 2
        # Add the intermediate node directionality constraint.
        constraint_intermediate_directionality(
            wm, nw, i, pipe_fr, pipe_to, des_pipe_fr, des_pipe_to, pump_fr, pump_to,
            regulator_fr, regulator_to, short_pipe_fr, short_pipe_to, valve_fr, valve_to)
    elseif length(reservoirs) > 0 && num_components == length(reservoirs)
        # Add the source node directionality constraint.
        constraint_source_directionality(
            wm, nw, i, pipe_fr, pipe_to, des_pipe_fr, des_pipe_to, pump_fr, pump_to,
            regulator_fr, regulator_to, short_pipe_fr, short_pipe_to, valve_fr, valve_to)
    elseif num_components == length(fixed_demands) && net_fixed_demand < 0.0
        # Add the source node directionality constraint.
        constraint_source_directionality(
            wm, nw, i, pipe_fr, pipe_to, des_pipe_fr, des_pipe_to, pump_fr, pump_to,
            regulator_fr, regulator_to, short_pipe_fr, short_pipe_to, valve_fr, valve_to)
    elseif num_components == length(fixed_demands) && net_fixed_demand > 0.0
        # Add the sink node directionality constraint.
        constraint_sink_directionality(
            wm, nw, i, pipe_fr, pipe_to, des_pipe_fr, des_pipe_to, pump_fr, pump_to,
            regulator_fr, regulator_to, short_pipe_fr, short_pipe_to, valve_fr, valve_to)
    end
end


### Tank Constraints ###
function constraint_tank_volume(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    if !(:time_step in keys(ref(wm, nw)))
        Memento.error(_LOGGER, "Tank states cannot be controlled without a time step.")
    end

    # Only set the tank state if the tank is nondispatchable.
    if !ref(wm, nw, :tank, i)["dispatchable"]
        tank = ref(wm, nw, :tank, i)
        initial_level = tank["init_level"]
        surface_area = 0.25 * pi * tank["diameter"]^2
        V_initial = surface_area * initial_level
        _initialize_con_dict(wm, :tank_volume, nw=nw)
        constraint_tank_volume_fixed(wm, nw, i, V_initial)
    end
end


function constraint_tank_volume(wm::AbstractWaterModel, i::Int, nw_1::Int, nw_2::Int)
    if !(:time_step in keys(ref(wm, nw_1)))
        Memento.error(_LOGGER, "Tank states cannot be controlled without a time step.")
    end

    # Only set the tank state if the tank is nondispatchable.
    if !ref(wm, nw_2, :tank, i)["dispatchable"]
        # TODO: What happens if a tank exists in nw_1 but not in nw_2? The index
        # "i" is assumed to be present in both when this constraint is applied.
        _initialize_con_dict(wm, :tank_volume, nw=nw_2)
        constraint_tank_volume(wm, nw_1, nw_2, i, ref(wm, nw_1, :time_step))
    end
end


### Pipe Constraints ###
function constraint_pipe_flow(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    _initialize_con_dict(wm, :pipe_flow, nw=nw, is_array=true)
    con(wm, nw, :pipe_flow)[a] = Array{JuMP.ConstraintRef}([])
    constraint_pipe_flow(wm, nw, a)
end


function constraint_pipe_head(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    _initialize_con_dict(wm, :pipe_head, nw=nw, is_array=true)
    con(wm, nw, :pipe_head)[a] = Array{JuMP.ConstraintRef}([])
    node_fr, node_to = ref(wm, nw, :pipe, a)["node_fr"], ref(wm, nw, :pipe, a)["node_to"]
    constraint_pipe_head(wm, nw, a, node_fr, node_to)
end


function constraint_pipe_head_loss(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    _initialize_con_dict(wm, :pipe_head_loss, nw=nw, is_array=true)
    con(wm, nw, :pipe_head_loss)[a] = Array{JuMP.ConstraintRef}([])
    node_fr, node_to = ref(wm, nw, :pipe, a)["node_fr"], ref(wm, nw, :pipe, a)["node_to"]
    exponent, L = ref(wm, nw, :alpha), ref(wm, nw, :pipe, a)["length"]
    r = minimum(ref(wm, nw, :resistance, a))
    constraint_pipe_head_loss(wm, nw, a, node_fr, node_to, exponent, L, r)
end


### Design Pipe Constraints ###
function constraint_on_off_pipe_flow_des(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    _initialize_con_dict(wm, :on_off_pipe_flow_des, nw=nw, is_array=true)
    con(wm, nw, :on_off_pipe_flow_des)[a] = Array{JuMP.ConstraintRef}([])
    resistances = ref(wm, nw, :resistance, a)
    constraint_on_off_pipe_flow_des(wm, nw, a, resistances)
end


function constraint_on_off_pipe_head_des(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    _initialize_con_dict(wm, :on_off_pipe_head_des, nw=nw, is_array=true)
    con(wm, nw, :on_off_pipe_head_des)[a] = Array{JuMP.ConstraintRef}([])
    node_fr = ref(wm, nw, :des_pipe, a)["node_fr"]
    node_to = ref(wm, nw, :des_pipe, a)["node_to"]
    constraint_on_off_pipe_head_des(wm, nw, a, node_fr, node_to)
end


function constraint_on_off_pipe_head_loss_des(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    _initialize_con_dict(wm, :on_off_pipe_head_loss_des, nw=nw, is_array=true)
    con(wm, nw, :on_off_pipe_head_loss_des)[a] = Array{JuMP.ConstraintRef}([])
    node_fr = ref(wm, nw, :des_pipe, a)["node_fr"]
    node_to = ref(wm, nw, :des_pipe, a)["node_to"]
    exponent, L = ref(wm, nw, :alpha), ref(wm, nw, :des_pipe, a)["length"]
    resist = ref(wm, nw, :resistance, a) # Dictionary of possibel resistances at `a`.
    constraint_on_off_pipe_head_loss_des(wm, nw, a, exponent, node_fr, node_to, L, resist)
end


### Pump Constraints ###
function constraint_on_off_pump_flow(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    node_fr, node_to = ref(wm, nw, :pump, a)["node_fr"], ref(wm, nw, :pump, a)["node_to"]
    _initialize_con_dict(wm, :on_off_pump_flow, nw=nw, is_array=true)
    con(wm, nw, :on_off_pump_flow)[a] = Array{JuMP.ConstraintRef}([])
    q_min_active = get(ref(wm, nw, :pump, a), "q_min_active", _q_eps)
    constraint_on_off_pump_flow(wm, nw, a, q_min_active)
end


function constraint_on_off_pump_head(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    node_fr, node_to = ref(wm, nw, :pump, a)["node_fr"], ref(wm, nw, :pump, a)["node_to"]
    _initialize_con_dict(wm, :on_off_pump_head, nw=nw, is_array=true)
    con(wm, nw, :on_off_pump_head)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_pump_head(wm, nw, a, node_fr, node_to)
end


function constraint_on_off_pump_head_gain(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    node_fr, node_to = ref(wm, nw, :pump, a)["node_fr"], ref(wm, nw, :pump, a)["node_to"]
    _initialize_con_dict(wm, :on_off_pump_head_gain, nw=nw, is_array=true)
    con(wm, nw, :on_off_pump_head_gain)[a] = Array{JuMP.ConstraintRef}([])
    coeffs = _get_function_from_head_curve(ref(wm, nw, :pump, a)["head_curve"])
    q_min_active = get(ref(wm, nw, :pump, a), "q_min_active", _q_eps)
    constraint_on_off_pump_head_gain(wm, nw, a, node_fr, node_to, coeffs, q_min_active)
end


### Short Pipe Constraints ###
function constraint_short_pipe_flow(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    _initialize_con_dict(wm, :short_pipe_flow, nw=nw, is_array=true)
    con(wm, nw, :short_pipe_flow)[a] = Array{JuMP.ConstraintRef}([])
    constraint_short_pipe_flow(wm, nw, a)
end


function constraint_short_pipe_head(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    node_fr, node_to = ref(wm, nw, :valve, a)["node_fr"], ref(wm, nw, :valve, a)["node_to"]
    _initialize_con_dict(wm, :short_pipe_head, nw=nw, is_array=true)
    con(wm, nw, :short_pipe_head)[a] = Array{JuMP.ConstraintRef}([])
    constraint_short_pipe_head(wm, nw, a, node_fr, node_to)
end


### Valve Constraints ###
function constraint_on_off_valve_flow(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    node_fr, node_to = ref(wm, nw, :valve, a)["node_fr"], ref(wm, nw, :valve, a)["node_to"]
    _initialize_con_dict(wm, :on_off_valve_flow, nw=nw, is_array=true)
    con(wm, nw, :on_off_valve_flow)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_valve_flow(wm, nw, a)
end


function constraint_on_off_valve_head(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    node_fr, node_to = ref(wm, nw, :valve, a)["node_fr"], ref(wm, nw, :valve, a)["node_to"]
    _initialize_con_dict(wm, :on_off_valve_head, nw=nw, is_array=true)
    con(wm, nw, :on_off_valve_head)[a] = Array{JuMP.ConstraintRef}([])
    constraint_on_off_valve_head(wm, nw, a, node_fr, node_to)
end


### Regulator Constraints ###
function constraint_on_off_regulator_flow(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    _initialize_con_dict(wm, :on_off_regulator_flow, nw=nw, is_array=true)
    con(wm, nw, :on_off_regulator_flow)[a] = Array{JuMP.ConstraintRef}([])
    q_min_active = get(ref(wm, nw, :regulator, a), "q_min_active", _q_eps)
    constraint_on_off_regulator_flow(wm, nw, a, q_min_active)
end


function constraint_on_off_regulator_head(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    _initialize_con_dict(wm, :on_off_regulator_head, nw=nw, is_array=true)
    con(wm, nw, :on_off_regulator_head)[a] = Array{JuMP.ConstraintRef}([])
    node_fr = ref(wm, nw, :regulator, a)["node_fr"]
    node_to = ref(wm, nw, :regulator, a)["node_to"]
    elevation = ref(wm, nw, :node, node_to)["elevation"]
    head_setting = elevation + ref(wm, nw, :regulator, a)["setting"]
    constraint_on_off_regulator_head(wm, nw, a, node_fr, node_to, head_setting)
end


### Pressure Reducing Valve Constraints ###
function constraint_prv_head_loss(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    node_fr = ref(wm, nw, :regulator, a)["node_fr"]
    node_to = ref(wm, nw, :regulator, a)["node_to"]

    # Add all common pressure reducing valve constraints.
    _initialize_con_dict(wm, :prv, nw=nw, is_array=true)
    con(wm, nw, :prv)[a] = Array{JuMP.ConstraintRef}([])
    constraint_prv_common(wm, nw, a, node_fr, node_to, h_prv)
end


### Pump Constraints ###
function constraint_pump_head_gain(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, force_on::Bool=false, kwargs...)
    # Get data common to all pump-related constraints.
    node_fr, node_to = ref(wm, nw, :pump, a)["node_fr"], ref(wm, nw, :pump, a)["node_to"]
    head_curve = ref(wm, nw, :pump, a)["head_curve"]
    coeffs = _get_function_from_head_curve(head_curve)

    # Add common constraints used to model pump behavior.
    _initialize_con_dict(wm, :pump, nw=nw, is_array=true)
    con(wm, nw, :pump)[a] = Array{JuMP.ConstraintRef}([])
    constraint_pump_common(wm, nw, a, node_fr, node_to, coeffs)

    # Add constraints that define head gain across the pump.
    _initialize_con_dict(wm, :head_gain, nw=nw, is_array=true)
    con(wm, nw, :head_gain)[a] = Array{JuMP.ConstraintRef}([])
    constraint_pump_head_gain(wm, nw, a, node_fr, node_to, coeffs)
    constraint_pump_head_gain_lb(wm, nw, a, node_fr, node_to, coeffs)
end


function constraint_head_difference_pump(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw)
    _initialize_con_dict(wm, :head_difference, nw=nw, is_array=true)
    con(wm, nw, :head_difference)[a] = Array{JuMP.ConstraintRef}([])
    node_fr, node_to = ref(wm, nw, :pump, a)["node_fr"], ref(wm, nw, :pump, a)["node_to"]
    constraint_head_difference_pump(wm, nw, a, node_fr, node_to)
end
