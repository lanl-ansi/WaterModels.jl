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

### Node Constraints ###
function constraint_source_head(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    h_src = ref(wm, nw, :reservoir, i)["elevation"]
    _initialize_con_dict(wm, :source_head, nw=nw)
    constraint_source_head(wm, nw, i, h_src)
end

function constraint_flow_conservation(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    arcs_from = ref(wm, nw, :node_arc_fr, i)
    arcs_to = ref(wm, nw, :node_arc_to, i)
    junctions = ref(wm, nw, :node_junction, i)
    reservoirs = ref(wm, nw, :node_reservoir, i)
    tanks = ref(wm, nw, :node_tank, i)
    demands = Dict{Int,Float64}(k => ref(wm, nw, :junction, k, "demand") for k in junctions)

    _initialize_con_dict(wm, :flow_conservation, nw=nw)
    constraint_flow_conservation(wm, nw, i, arcs_from, arcs_to, reservoirs, tanks, demands)
end

function constraint_sink_flow(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    _initialize_con_dict(wm, :sink_flow, nw=nw)
    constraint_sink_flow(wm, nw, i, ref(wm, nw, :link))
end

function constraint_source_flow(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    _initialize_con_dict(wm, :source_flow, nw=nw)
    constraint_source_flow(wm, nw, i, ref(wm, nw, :link))
end

### Tank Constraints ###
function constraint_link_volume(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    tank = ref(wm, nw, :tank, i)
    elevation = ref(wm, nw, :node, tank["tank_node"])["elevation"]
    surface_area = 0.25 * pi * tank["diameter"]^2

    _initialize_con_dict(wm, :link_volume, nw=nw)
    constraint_link_volume(wm, nw, i, elevation, surface_area)
end

function constraint_tank_state(wm::AbstractWaterModel, i::Int; nw::Int=wm.cnw)
    if "hydraulic_timestep" in keys(ref(wm, nw, :option, "time"))
        time_step_int = ref(wm, nw, :option, "time")["hydraulic_timestep"]
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
    if "hydraulic_timestep" in keys(ref(wm, nw_2, :option, "time"))
        time_step_int = ref(wm, nw_2, :option, "time")["hydraulic_timestep"]
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
function constraint_head_loss_pipe(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    alpha = ref(wm, nw, :alpha)
    pipe = ref(wm, nw, :pipe, a)
    r = minimum(ref(wm, nw, :resistance, a))

    _initialize_con_dict(wm, :head_loss, nw=nw, is_array=true)
    con(wm, nw, :head_loss)[a] = Array{JuMP.ConstraintRef}([])
    constraint_flow_direction_selection(wm, a; nw=nw, kwargs...)
    constraint_head_difference(wm, a; nw=nw, kwargs...)
    constraint_head_loss_ub_pipe(wm, a; nw=nw, kwargs...)
    constraint_head_loss_pipe(wm, nw, a, alpha, pipe["node_fr"], pipe["node_to"], pipe["length"], r)
end

function constraint_head_difference(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw)
    node_fr = ref(wm, nw, :link, a)["node_fr"]
    head_fr = nothing

    for rid in ref(wm, nw, :node_reservoir, node_fr)
        # TODO: This is a good place to check these are consistent.
        head_fr = ref(wm, nw, :reservoir, rid)["elevation"]
    end

    node_to = ref(wm, nw, :link, a)["node_to"]
    head_to = nothing

    for rid in ref(wm, nw, :node_reservoir, node_to)
        # TODO: This is a good place to check these are consistent
        head_to = ref(wm, nw, :reservoir, rid)["elevation"]
    end

    constraint_head_difference(wm, nw, a, node_fr, node_to, head_fr, head_to)
end

function constraint_flow_direction_selection(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw)
    constraint_flow_direction_selection(wm, nw, a)
end

function constraint_head_loss_ub_pipe(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw)
    alpha = ref(wm, nw, :alpha)
    L = ref(wm, nw, :pipe, a)["length"]
    r = maximum(ref(wm, nw, :resistance, a))
    constraint_head_loss_ub_pipe(wm, nw, a, alpha, L, r)
end

function constraint_head_loss_pipe_ne(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    alpha = ref(wm, nw, :alpha)
    pipe = ref(wm, nw, :pipe, a)
    resistances = ref(wm, nw, :resistance, a)

    _initialize_con_dict(wm, :head_loss, nw=nw, is_array=true)
    con(wm, nw, :head_loss)[a] = Array{JuMP.ConstraintRef}([])
    constraint_flow_direction_selection_ne(wm, a; nw=nw, kwargs...)
    constraint_head_difference(wm, a; nw=nw, kwargs...)
    constraint_head_loss_pipe_ne(wm, nw, a, alpha, pipe["node_fr"], pipe["node_to"], pipe["length"], resistances)
    constraint_head_loss_ub_pipe_ne(wm, a; nw=nw, kwargs...)
    constraint_resistance_selection_ne(wm, a; nw=nw, kwargs...)
end

function constraint_resistance_selection_ne(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    pipe_resistances = ref(wm, nw, :resistance, a)
    constraint_resistance_selection_ne(wm, nw, a, pipe_resistances; kwargs...)
end

function constraint_flow_direction_selection_ne(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw)
    pipe_resistances = ref(wm, nw, :resistance, a)
    constraint_flow_direction_selection_ne(wm, nw, a, pipe_resistances)
end

function constraint_head_loss_ub_pipe_ne(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw)
    alpha = ref(wm, nw, :alpha)
    L = ref(wm, nw, :pipe, a)["length"]
    pipe_resistances = ref(wm, nw, :resistance, a)
    constraint_head_loss_ub_pipe_ne(wm, nw, a, alpha, L, pipe_resistances)
end

### Check Valve Constraints ###
function constraint_check_valve(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw)
    node_fr = ref(wm, nw, :link, a)["node_fr"]
    node_to = ref(wm, nw, :link, a)["node_to"]

    _initialize_con_dict(wm, :check_valve, nw=nw, is_array=true)
    con(wm, nw, :check_valve)[a] = Array{JuMP.ConstraintRef}([])
    constraint_check_valve(wm, nw, a, node_fr, node_to)
end

function constraint_head_loss_check_valve(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw)
    node_fr = ref(wm, nw, :link, a)["node_fr"]
    node_to = ref(wm, nw, :link, a)["node_to"]
    L = ref(wm, nw, :pipe, a)["length"]
    r = minimum(ref(wm, nw, :resistance, a))

    _initialize_con_dict(wm, :head_loss, nw=nw, is_array=true)
    con(wm, nw, :head_loss)[a] = Array{JuMP.ConstraintRef}([])
    constraint_head_loss_check_valve(wm, nw, a, node_fr, node_to, L, r)
end

### Pump Constraints ###
function constraint_head_gain_pump(wm::AbstractWaterModel, a::Int; nw::Int=wm.cnw, force_on::Bool=false, kwargs...)
    node_fr = ref(wm, nw, :pump, a)["node_fr"]
    node_to = ref(wm, nw, :pump, a)["node_to"]
    pump_curve = ref(wm, nw, :pump, a)["pump_curve"]
    curve_fun = _get_function_from_pump_curve(pump_curve)

    _initialize_con_dict(wm, :head_gain, nw=nw, is_array=true)
    con(wm, nw, :head_gain)[a] = Array{JuMP.ConstraintRef}([])
    force_on ? constraint_head_gain_pump_on(wm, nw, a, node_fr, node_to, curve_fun) :
        constraint_head_gain_pump(wm, nw, a, node_fr, node_to, curve_fun)
end

function constraint_pump_control(wm::AbstractWaterModel, a::Int, nw_1::Int, nw_2::Int)
    pump = ref(wm, nw_2, :pump, a)

    if "control" in keys(pump)
        comps = Array{Pair{String, Int64}, 1}()
        lt, gt = [nothing, nothing]

        for (control_name, control) in pump["control"]
            action = control["action"]
            condition = control["condition"]
            comps = vcat((condition["node_type"], condition["node_id"]), comps)

            if condition["operator"] == "<="
                lt = condition["threshold"]
            elseif condition["operator"] == ">="
                gt = condition["threshold"]
            end
        end

        if all(y->y==comps[1], comps) && lt <= gt && comps[1][1] == "tank"
            elevation = ref(wm, nw_2, :node, comps[1][2])["elevation"]
            _initialize_con_dict(wm, :pump_control, nw=nw_2)
            constraint_pump_control_tank(wm, nw_1, nw_2, a, comps[1][2], lt, gt, elevation)
        else
            Memento.error(_LOGGER, "Can't handle control condition.")
        end
    end
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
