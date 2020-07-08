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

### Nodal Constraints ###
function constraint_flow_conservation(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    links_from = ref(wm, nw, :node_link_fr, i)
    links_to = ref(wm, nw, :node_link_to, i)
    junctions = ref(wm, nw, :node_junction, i)
    reservoirs = ref(wm, nw, :node_reservoir, i)
    tanks = ref(wm, nw, :node_tank, i)
    demands = Dict{Int,Float64}(k => ref(wm, nw, :junction, k, "demand") for k in junctions)

    _initialize_con_dict(wm, :flow_conservation, nw=nw)
    constraint_flow_conservation(wm, nw, i, links_from, links_to, reservoirs, tanks, demands)
end

### Junction Constraints ###
function constraint_sink_flow(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    links_from = ref(wm, nw, :node_link_fr, i)
    links_to = ref(wm, nw, :node_link_to, i)
    _initialize_con_dict(wm, :sink_flow, nw=nw)
    constraint_sink_flow(wm, nw, i, links_from, links_to)
end

### Reservoir Constraints ###
function constraint_source_head(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    h_src = ref(wm, nw, :reservoir, i)["elevation"]
    _initialize_con_dict(wm, :source_head, nw=nw)
    constraint_source_head(wm, nw, i, h_src)
end

function constraint_source_flow(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    links_from = ref(wm, nw, :node_link_fr, i)
    links_to = ref(wm, nw, :node_link_to, i)
    _initialize_con_dict(wm, :source_flow, nw=nw)
    constraint_source_flow(wm, nw, i, links_from, links_to)
end

### Tank Constraints ###
function constraint_volume(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    tank = ref(wm, nw, :tank, i)
    elevation = ref(wm, nw, :node, tank["tank_node"])["elevation"]
    surface_area = 0.25 * pi * tank["diameter"]^2

    _initialize_con_dict(wm, :volume, nw=nw)
    constraint_volume(wm, nw, i, elevation, surface_area)
end

function constraint_tank_state(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    if "hydraulic_timestep" in keys(wm.ref[:option]["time"])
        time_step_int = wm.ref[:option]["time"]["hydraulic_timestep"]
        time_step = convert(Float64, time_step_int)
    else
        Memento.error(_LOGGER, "Tank states cannot be controlled without a time step.")
    end

    tank = ref(wm, nw, :tank, i)
    initial_level = tank["init_level"]
    surface_area = 0.25 * pi * tank["diameter"]^2
    V_initial = surface_area * initial_level

    _initialize_con_dict(wm, :tank_state, nw=nw)
    constraint_tank_state_initial(wm, nw, i, V_initial, time_step)
end

function constraint_tank_state(wm::AbstractWaterModel, i::Int, nw_1::Int, nw_2::Int)
    if "hydraulic_timestep" in keys(wm.ref[:option]["time"])
        time_step_int = wm.ref[:option]["time"]["hydraulic_timestep"]
        time_step = convert(Float64, time_step_int)
    else
        Memento.error(_LOGGER, "Tank states cannot be controlled without a time step.")
    end

    # TODO: What happens if a tank exists in nw_1 but not in nw_2? The index
    # "i" is assumed to be present in both when this constraint is applied.
    _initialize_con_dict(wm, :tank_state, nw=nw_2)
    constraint_tank_state(wm, nw_1, nw_2, i, time_step)
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
    node_fr = ref(wm, nw, :shutoff_valve, a)["node_fr"]
    node_to = ref(wm, nw, :shutoff_valve, a)["node_to"]
    alpha, L = ref(wm, nw, :alpha), ref(wm, nw, :shutoff_valve, a)["length"]
    r = maximum(ref(wm, nw, :resistance, a))

    ## Since shutoff valves exist along pipes, add all common pipe constraints.
    #_initialize_con_dict(wm, :pipe, nw=nw, is_array=true)
    #con(wm, nw, :pipe)[a] = Array{JuMP.ConstraintRef}([])
    #constraint_pipe_common(wm, nw, a, node_fr, node_to, alpha, L, r)

    # Add all common shutoff valve constraints.
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
    h_lb = ref(wm, nw, :node, node_to)["elevation"]
    h_prv = h_lb + ref(wm, nw, :pressure_reducing_valve, a)["setting"]

    # Add all common pressure reducing valve constraints.
    _initialize_con_dict(wm, :prv, nw=nw, is_array=true)
    con(wm, nw, :prv)[a] = Array{JuMP.ConstraintRef}([])
    constraint_prv_common(wm, nw, a, node_fr, node_to, h_prv)
end

### Pump Constraints ###
function constraint_pump_head_gain(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, force_on::Bool=false, kwargs...)
    # Get data common to all pump-related constraints.
    node_fr, node_to = ref(wm, nw, :pump, a)["node_fr"], ref(wm, nw, :pump, a)["node_to"]
    pump_curve = ref(wm, nw, :pump, a)["pump_curve"]
    coeffs = _get_function_from_pump_curve(pump_curve)

    # Add common constraints used to model pump behavior.
    _initialize_con_dict(wm, :pump, nw=nw, is_array=true)
    con(wm, nw, :pump)[a] = Array{JuMP.ConstraintRef}([])
    constraint_pump_common(wm, nw, a, node_fr, node_to, coeffs)

    # Add constraints that define head gain across the pump.
    _initialize_con_dict(wm, :head_gain, nw=nw, is_array=true)
    con(wm, nw, :head_gain)[a] = Array{JuMP.ConstraintRef}([])
    constraint_pump_head_gain(wm, nw, a, node_fr, node_to, coeffs)
end

function constraint_head_difference_pump(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw)
    _initialize_con_dict(wm, :head_difference, nw=nw, is_array=true)
    con(wm, nw, :head_difference)[a] = Array{JuMP.ConstraintRef}([])
    node_fr, node_to = ref(wm, nw, :pump, a)["node_fr"], ref(wm, nw, :pump, a)["node_to"]
    constraint_head_difference_pump(wm, nw, a, node_fr, node_to)
end

function constraint_pump_control(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw)
    pump = ref(wm, nw, :pump, a)

    if "initial_status" in keys(pump)
        initial_status = uppercase(pump["initial_status"]) != "CLOSED"
        _initialize_con_dict(wm, :pump_control, nw=nw)
        constraint_pump_control_initial(wm, nw, a, initial_status)
    else
        Memento.error(_LOGGER, "Initial status not found for pump.")
    end
end
