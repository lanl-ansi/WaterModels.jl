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
    check_valve_fr = _collect_comps_fr(wm, i, :check_valve; nw=nw)
    check_valve_to = _collect_comps_to(wm, i, :check_valve; nw=nw)
    pipe_fr = _collect_comps_fr(wm, i, :pipe; nw=nw)
    pipe_to = _collect_comps_to(wm, i, :pipe; nw=nw)
    pump_fr = _collect_comps_fr(wm, i, :pump; nw=nw)
    pump_to = _collect_comps_to(wm, i, :pump; nw=nw)
    pressure_reducing_valve_fr = _collect_comps_fr(wm, i, :pressure_reducing_valve; nw=nw)
    pressure_reducing_valve_to = _collect_comps_to(wm, i, :pressure_reducing_valve; nw=nw)
    shutoff_valve_fr = _collect_comps_fr(wm, i, :shutoff_valve; nw=nw)
    shutoff_valve_to = _collect_comps_to(wm, i, :shutoff_valve; nw=nw)

    # Collect various indices for node-type components connected to node `i`.
    reservoirs = ref(wm, nw, :node_reservoir, i) # Reservoirs attached to node `i`.
    tanks = ref(wm, nw, :node_tank, i) # Tanks attached to node `i`.
    junctions = ref(wm, nw, :node_junction, i) # Junctions attached to node `i`.

    # Sum the constant demands required at node `i`.
    nondispatchable_junctions = filter(j -> j in ids(wm, nw, :nondispatchable_junction), junctions)
    fixed_demands = [ref(wm, nw, :nondispatchable_junction, j)["demand"] for j in nondispatchable_junctions]
    net_fixed_demand = length(fixed_demands) > 0 ? sum(fixed_demands) : 0.0

    # Get the indices of dispatchable junctions connected to node `i`.
    dispatchable_junctions = filter(j -> j in ids(wm, nw, :dispatchable_junction), junctions)

    # Initialize the flow conservation constraint dictionary entry.
    _initialize_con_dict(wm, :flow_conservation, nw=nw)

    # Add the flow conservation constraint.
    constraint_flow_conservation(
        wm, nw, i, check_valve_fr, check_valve_to, pipe_fr, pipe_to, pump_fr, pump_to,
        pressure_reducing_valve_fr, pressure_reducing_valve_to, shutoff_valve_fr,
        shutoff_valve_to, reservoirs, tanks, dispatchable_junctions, net_fixed_demand)
end


function constraint_node_directionality(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    # Collect various indices for edge-type components connected to node `i`.
    check_valve_fr = _collect_comps_fr(wm, i, :check_valve; nw=nw)
    check_valve_to = _collect_comps_to(wm, i, :check_valve; nw=nw)
    pipe_fr = _collect_comps_fr(wm, i, :pipe; nw=nw)
    pipe_to = _collect_comps_to(wm, i, :pipe; nw=nw)
    pump_fr = _collect_comps_fr(wm, i, :pump; nw=nw)
    pump_to = _collect_comps_to(wm, i, :pump; nw=nw)
    pressure_reducing_valve_fr = _collect_comps_fr(wm, i, :pressure_reducing_valve; nw=nw)
    pressure_reducing_valve_to = _collect_comps_to(wm, i, :pressure_reducing_valve; nw=nw)
    shutoff_valve_fr = _collect_comps_fr(wm, i, :shutoff_valve; nw=nw)
    shutoff_valve_to = _collect_comps_to(wm, i, :shutoff_valve; nw=nw)

    # Get the number of nodal components attached to node `i`.
    junctions = ref(wm, nw, :node_junction)
    tanks = ref(wm, nw, :node_tank)
    reservoirs = ref(wm, nw, :node_reservoir)
    num_components = length(junctions) + length(tanks) + length(reservoirs)

    # Get the in degree of node `i`.
    in_length = length(check_valve_to) + length(pipe_to) + length(pump_to) +
                length(pressure_reducing_valve_to) + length(shutoff_valve_to)

    # Get the out degree of node `i`.
    out_length = length(check_valve_fr) + length(pipe_fr) + length(pump_fr) +
                 length(pressure_reducing_valve_fr) + length(shutoff_valve_fr)

    # Check if node directionality constraints should be added.
    if num_components == 0 && in_length + out_length == 2
        # Initialize the node directionality constraint dictionary entry.
        _initialize_con_dict(wm, :node_directionality, nw=nw)

        # Add the node directionality constraint.
        constraint_node_directionality(
            wm, nw, i, check_valve_fr, check_valve_to, pipe_fr, pipe_to, pump_fr, pump_to,
            pressure_reducing_valve_fr, pressure_reducing_valve_to, shutoff_valve_fr,
            shutoff_valve_to)
    end
end


### Junction Constraints ###
function constraint_sink_directionality(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    # Collect various indices for edge-type components connected to node `i`.
    check_valve_fr = _collect_comps_fr(wm, i, :check_valve; nw=nw)
    check_valve_to = _collect_comps_to(wm, i, :check_valve; nw=nw)
    pipe_fr = _collect_comps_fr(wm, i, :pipe; nw=nw)
    pipe_to = _collect_comps_to(wm, i, :pipe; nw=nw)
    pump_fr = _collect_comps_fr(wm, i, :pump; nw=nw)
    pump_to = _collect_comps_to(wm, i, :pump; nw=nw)
    pressure_reducing_valve_fr = _collect_comps_fr(wm, i, :pressure_reducing_valve; nw=nw)
    pressure_reducing_valve_to = _collect_comps_to(wm, i, :pressure_reducing_valve; nw=nw)
    shutoff_valve_fr = _collect_comps_fr(wm, i, :shutoff_valve; nw=nw)
    shutoff_valve_to = _collect_comps_to(wm, i, :shutoff_valve; nw=nw)

    # Initialize the sink flow constraint dictionary entry.
    _initialize_con_dict(wm, :sink_directionality, nw=nw)

    # Add the sink flow directionality constraint.
    constraint_sink_directionality(
        wm, nw, i, check_valve_fr, check_valve_to, pipe_fr, pipe_to, pump_fr, pump_to,
        pressure_reducing_valve_fr, pressure_reducing_valve_to, shutoff_valve_fr,
        shutoff_valve_to)
end


### Reservoir Constraints ###
function constraint_reservoir_head(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    # Only fix the reservoir head if the reservoir is nondispatchable.
    if !ref(wm, nw, :reservoir, i)["dispatchable"]
        node_index = ref(wm, nw, :reservoir, i)["node"]
        head = ref(wm, nw, :node, node_index)["head"]
        _initialize_con_dict(wm, :reservoir_head, nw=nw)
        constraint_reservoir_head(wm, nw, node_index, head)
    end
end


function constraint_source_directionality(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    # Collect various indices for edge-type components connected to node `i`.
    check_valve_fr = _collect_comps_fr(wm, i, :check_valve; nw=nw)
    check_valve_to = _collect_comps_to(wm, i, :check_valve; nw=nw)
    pipe_fr = _collect_comps_fr(wm, i, :pipe; nw=nw)
    pipe_to = _collect_comps_to(wm, i, :pipe; nw=nw)
    pump_fr = _collect_comps_fr(wm, i, :pump; nw=nw)
    pump_to = _collect_comps_to(wm, i, :pump; nw=nw)
    pressure_reducing_valve_fr = _collect_comps_fr(wm, i, :pressure_reducing_valve; nw=nw)
    pressure_reducing_valve_to = _collect_comps_to(wm, i, :pressure_reducing_valve; nw=nw)
    shutoff_valve_fr = _collect_comps_fr(wm, i, :shutoff_valve; nw=nw)
    shutoff_valve_to = _collect_comps_to(wm, i, :shutoff_valve; nw=nw)

    # Initialize the source flow constraint dictionary entry.
    _initialize_con_dict(wm, :source_directionality, nw=nw)

    # Add the source flow directionality constraint.
    constraint_source_directionality(
        wm, nw, i, check_valve_fr, check_valve_to, pipe_fr, pipe_to, pump_fr,
        pump_to, pressure_reducing_valve_fr, pressure_reducing_valve_to, shutoff_valve_fr,
        shutoff_valve_to)
end


### Tank Constraints ###
function constraint_tank_state(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    if !(:time_step in keys(ref(wm, nw)))
        Memento.error(_LOGGER, "Tank states cannot be controlled without a time step.")
    end

    # Only set the tank state if the tank is nondispatchable.
    if !ref(wm, nw, :tank, i)["dispatchable"]
        tank = ref(wm, nw, :tank, i)
        initial_level = tank["init_level"]
        surface_area = 0.25 * pi * tank["diameter"]^2
        V_initial = surface_area * initial_level
        _initialize_con_dict(wm, :tank_state, nw=nw)
        constraint_tank_state_initial(wm, nw, i, V_initial)
    end
end


function constraint_tank_state(wm::AbstractWaterModel, i::Int, nw_1::Int, nw_2::Int)
    if !(:time_step in keys(ref(wm, nw_1)))
        Memento.error(_LOGGER, "Tank states cannot be controlled without a time step.")
    end

    # Only set the tank state if the tank is nondispatchable.
    if !ref(wm, nw_1, :tank, i)["dispatchable"]
        # TODO: What happens if a tank exists in nw_1 but not in nw_2? The index
        # "i" is assumed to be present in both when this constraint is applied.
        _initialize_con_dict(wm, :tank_state, nw=nw_2)
        constraint_tank_state(wm, nw_1, nw_2, i, ref(wm, nw_1, :time_step))
    end
end


### Pipe Constraints ###
function constraint_pipe_head_loss(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    node_fr, node_to = ref(wm, nw, :pipe, a)["node_fr"], ref(wm, nw, :pipe, a)["node_to"]
    alpha, L = ref(wm, nw, :alpha), ref(wm, nw, :pipe, a)["length"]
    r = minimum(ref(wm, nw, :resistance, a))

    # Add common constraints used to model pump behavior.
    _initialize_con_dict(wm, :pipe, nw=nw, is_array=true)
    con(wm, nw, :pipe)[a] = Array{JuMP.ConstraintRef}([])
    constraint_pipe_common(wm, nw, a, node_fr, node_to, alpha, L, r)

    # Add constraints related to modeling head loss along pipe.
    _initialize_con_dict(wm, :head_loss, nw=nw, is_array=true)
    con(wm, nw, :head_loss)[a] = Array{JuMP.ConstraintRef}([])
    constraint_pipe_head_loss(wm, nw, a, node_fr, node_to, alpha, L, r)
    constraint_pipe_head_loss_ub(wm, nw, a, alpha, L, r)
end


function constraint_pipe_head_loss_des(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    alpha = ref(wm, nw, :alpha)
    pipe = ref(wm, nw, :pipe, a)
    resistances = ref(wm, nw, :resistance, a)
    i, j = pipe["node_fr"], pipe["node_to"]

    _initialize_con_dict(wm, :head_loss, nw=nw, is_array=true)
    con(wm, nw, :head_loss)[a] = Array{JuMP.ConstraintRef}([])

    constraint_flow_direction_selection_des(wm, a; nw=nw, kwargs...)
    constraint_pipe_head_loss_des(wm, nw, a, alpha, i, j, pipe["length"], resistances)
    constraint_pipe_head_loss_ub_des(wm, a; nw=nw, kwargs...)
    constraint_resistance_selection_des(wm, a; nw=nw, kwargs...)
end


function constraint_resistance_selection_des(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    pipe_resistances = ref(wm, nw, :resistance, a)
    constraint_resistance_selection_des(wm, nw, a, pipe_resistances; kwargs...)
end


function constraint_flow_direction_selection_des(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw)
    pipe_resistances = ref(wm, nw, :resistance, a)
    constraint_flow_direction_selection_des(wm, nw, a, pipe_resistances)
end


function constraint_pipe_head_loss_ub_des(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw)
    alpha, L = [ref(wm, nw, :alpha), ref(wm, nw, :pipe, a)["length"]]
    pipe_resistances = ref(wm, nw, :resistance, a)
    constraint_pipe_head_loss_ub_des(wm, nw, a, alpha, L, pipe_resistances)
end


### Check Valve Constraints ###
function constraint_check_valve_head_loss(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    node_fr = ref(wm, nw, :check_valve, a)["node_fr"]
    node_to = ref(wm, nw, :check_valve, a)["node_to"]
    alpha, L = ref(wm, nw, :alpha), ref(wm, nw, :check_valve, a)["length"]
    r = maximum(ref(wm, nw, :resistance, a))

    # Add all common check valve constraints.
    _initialize_con_dict(wm, :check_valve, nw=nw, is_array=true)
    con(wm, nw, :check_valve)[a] = Array{JuMP.ConstraintRef}([])
    constraint_check_valve_common(wm, nw, a, node_fr, node_to)

    # Add constraints that describe head loss relationships.
    _initialize_con_dict(wm, :head_loss, nw=nw, is_array=true)
    alpha, L = [ref(wm, nw, :alpha), ref(wm, nw, :check_valve, a)["length"]]
    r = minimum(ref(wm, nw, :resistance, a))
    con(wm, nw, :head_loss)[a] = Array{JuMP.ConstraintRef}([])
    constraint_check_valve_head_loss(wm, nw, a, node_fr, node_to, L, r)
    constraint_head_loss_ub_cv(wm, nw, a, alpha, L, r)
end


### Shutoff Valve Constraints ###
function constraint_shutoff_valve_head_loss(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    # Add all common shutoff valve constraints.
    node_fr = ref(wm, nw, :shutoff_valve, a)["node_fr"]
    node_to = ref(wm, nw, :shutoff_valve, a)["node_to"]
    _initialize_con_dict(wm, :sv, nw=nw, is_array=true)
    con(wm, nw, :sv)[a] = Array{JuMP.ConstraintRef}([])
    constraint_sv_common(wm, nw, a, node_fr, node_to)

    # Add constraints that describe head loss relationships.
    _initialize_con_dict(wm, :head_loss, nw=nw, is_array=true)
    alpha, L = ref(wm, nw, :alpha), ref(wm, nw, :shutoff_valve, a)["length"]
    r = minimum(ref(wm, nw, :resistance, a))
    con(wm, nw, :head_loss)[a] = Array{JuMP.ConstraintRef}([])
    constraint_shutoff_valve_head_loss(wm, nw, a, node_fr, node_to, L, r)
    constraint_shutoff_valve_head_loss_ub(wm, nw, a, alpha, L, r)
end


### Pressure Reducing Valve Constraints ###
function constraint_prv_head_loss(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    node_fr = ref(wm, nw, :pressure_reducing_valve, a)["node_fr"]
    node_to = ref(wm, nw, :pressure_reducing_valve, a)["node_to"]
    elev = ref(wm, nw, :node, node_to)["elevation"]
    h_prv = elev + ref(wm, nw, :pressure_reducing_valve, a)["setting"]

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
