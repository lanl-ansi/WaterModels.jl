# Constraint templates help simplify data wrangling across multiple Water Flow
# formulations by providing an abstraction layer between the network data and
# network constraint definitions. The constraint template's job is to extract
# the required parameters from a given network data structure and pass the data
# as named arguments to the Water Flow formulations.
#
# Constraint templates should always be defined over "GenericWaterModel" and
# should never refer to model variables.


### Node Constraints ###
function constraint_flow_conservation(wm::GenericWaterModel, i::Int; nw::Int=wm.cnw)
    # Create the constraint dictionary if necessary.
    if !haskey(con(wm, nw), :flow_conservation)
        con(wm, nw)[:flow_conservation] = Dict{Int, JuMP.ConstraintRef}()
    end

    node_junctions = ref(wm, nw, :node_junctions, i)
    node_arcs_fr = ref(wm, nw, :node_arcs_fr, i)
    node_arcs_to = ref(wm, nw, :node_arcs_to, i)
    node_reservoirs = ref(wm, nw, :node_reservoirs, i)
    node_tanks = ref(wm, nw, :node_tanks, i)

    node_demands = Dict(k => ref(wm, nw, :junctions, k, "demand") for k in node_junctions)

    constraint_flow_conservation(wm, nw, i, node_arcs_fr, node_arcs_to, node_reservoirs, node_tanks, node_demands)
end


function constraint_sink_flow(wm::GenericWaterModel, i::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, n), :directed_sink_flow)
        con(wm, n)[:directed_sink_flow] = Dict{Int, JuMP.ConstraintRef}()
    end

    links = ref(wm, nw, :links)
    constraint_sink_flow(wm, nw, i, links)
end

function constraint_source_flow(wm::GenericWaterModel, i::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, n), :directed_source_flow)
        con(wm, n)[:directed_source_flow] = Dict{Int, JuMP.ConstraintRef}()
    end

    links = ref(wm, nw, :links)
    constraint_sink_flow(wm, nw, i, links)
end



### Junction Constraints ###



### Reservoir Constraints ###



### Tank Constraints ###
function constraint_link_volume(wm::GenericWaterModel, i::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, nw), :link_volume)
        con(wm, nw)[:link_volume] = Dict{Int, JuMP.ConstraintRef}()
    end

    tank = ref(wm, nw, :tanks, i)
    elevation = ref(wm, nw, :nodes, tank["tank_node"])["elevation"]
    surface_area = 0.25 * pi * tank["diameter"]^2

    constraint_link_volume(wm, nw, i, elevation, surface_area)
end

""
function constraint_pump_control(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, nw), :pump_control)
        con(wm, nw)[:pump_control] = Dict{Int, JuMP.ConstraintRef}()
    end

    pump = ref(wm, nw, :pumps, a)

    if "controls" in keys(pump)
        comps = Array{Pair{String, Int64}, 1}()
        lt, gt = [nothing, nothing]

        for (control_name, control) in pump["controls"]
            action = control["action"]
            condition = control["condition"]
            comps = vcat((condition["node_type"], condition["node_id"]), comps)

            if condition["operator"] == "<="
                lt = condition["threshold"]
            elseif condition["operator"] == ">="
                gt = condition["threshold"]
            end
        end

        same_comp = all(y->y==comps[1], comps)

        if same_comp && lt <= gt && comps[1][1] == "tanks"
            elevation = ref(wm, nw, :nodes, comps[1][2])["elevation"]
            constraint_pump_control_tank(wm, nw, a, comps[1][2], lt, gt, elevation)
        else
            Memento.error(_LOGGER, "Can't handle control condition.")
        end
    end
end

""
function constraint_tank_state(wm::GenericWaterModel, i::Int; nw::Int=wm.cnw)
    tank = ref(wm, nw, :tanks, i)
    time_step = ref(wm, nw, :options, "time")["hydraulic_timestep"]

    if time_step <= 0.0
        Memento.error(_LOGGER, "Tank states cannot be controlled without a time step.")
    end

    initial_level = tank["init_level"]
    surface_area = 0.25 * pi * tank["diameter"]^2
    V_initial = surface_area * initial_level

    constraint_tank_state_initial(wm, nw, i, V_initial, convert(Float64, time_step))
end

function constraint_tank_state(wm::GenericWaterModel, i::Int, nw_1::Int, nw_2::Int)
    tank = ref(wm, nw_2, :tanks, i)

    if haskey(ref(wm, nw_2, :options)["time"], "hydraulic_timestep")
        time_step = ref(wm, nw_2, :options)["time"]["hydraulic_timestep"]
    else
        Memento.error(_LOGGER, "Tank states cannot be controlled outside of a time series.")
    end

    constraint_tank_state(wm, nw_1, nw_2, i, convert(Float64, time_step))
end



### Link Constraints ###
function constraint_link_flow(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, nw), :link_directed_flow)
        con(wm, nw)[:link_directed_flow] = Dict{Int, JuMP.ConstraintRef}()
    end

    constraint_link_flow(wm, nw, a)
end

function constraint_link_flow_ne(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    constraint_link_flow_ne(wm, nw, a)
end


### Pipe Constraints ###
function constraint_potential_loss_pipe(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    alpha = ref(wm, nw, :alpha)

    pipe = ref(wm, nw, :pipes, a)
    r_min = minimum(ref(wm, nw, :resistance, a))

    constraint_potential_loss_pipe(wm, nw, a, alpha, pipe["f_id"], pipe["t_id"], pipe["length"], r_min)

    constraint_head_difference(wm, a; nw=nw, kwargs...)
    constraint_flow_direction_selection(wm, a; nw=nw, kwargs...)
    constraint_potential_loss_ub_pipe(wm, a; nw=nw, kwargs...)
end

function constraint_head_difference(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, nw), :head_difference_1)
        con(wm, nw)[:head_difference_1] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, nw)[:head_difference_2] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, nw)[:head_difference_3] = Dict{Int, JuMP.ConstraintRef}()
    end

    f_id = ref(wm, nw, :links, a)["f_id"]

    head_fr = nothing
    for rid in ref(wm, nw, :node_reservoirs, f_id)
        #TODO this is a good place to check these are consistent
        head_fr = ref(wm, nw, :reservoirs, rid)["head"]
    end

    t_id = ref(wm, nw, :links, a)["t_id"]

    head_to = nothing
    for rid in ref(wm, nw, :node_reservoirs, t_id)
        #TODO this is a good place to check these are consistent
        head_to = ref(wm, nw, :reservoirs, rid)["head"]
    end

    constraint_head_difference(wm, nw, a, f_id, t_id, head_fr, head_to)
end

function constraint_flow_direction_selection(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, nw), :flow_direction_selection_n)
        con(wm, nw)[:flow_direction_selection_n] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, nw)[:flow_direction_selection_p] = Dict{Int, JuMP.ConstraintRef}()
    end

    constraint_flow_direction_selection(wm, nw, a)
end

function constraint_potential_loss_ub_pipe(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, nw), :directed_potential_loss_ub_pipe_n)
        con(wm, nw)[:directed_potential_loss_ub_pipe_p] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, nw)[:directed_potential_loss_ub_pipe_n] = Dict{Int, JuMP.ConstraintRef}()
    end

    alpha = ref(wm, nw, :alpha)
    len = ref(wm, nw, :pipes, a)["length"]
    r_max = maximum(ref(wm, nw, :resistance, a))

    constraint_potential_loss_ub_pipe(wm, nw, a, alpha, len, r_max)
end


function constraint_potential_loss_pipe_ne(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    alpha = ref(wm, nw, :alpha)

    pipe = ref(wm, nw, :pipes, a)
    pipe_resistances = ref(wm, nw, :resistance, a)

    constraint_potential_loss_pipe_ne(wm, nw, a, alpha, pipe["f_id"], pipe["t_id"], pipe["length"], pipe_resistances)

    constraint_head_difference(wm, a; nw=nw, kwargs...)
    constraint_flow_direction_selection_ne(wm, a; nw=nw, kwargs...)
    constraint_potential_loss_ub_pipe_ne(wm, a; nw=nw, kwargs...)
end


function constraint_flow_direction_selection_ne(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, nw), :flow_direction_selection_ne_n)
        con(wm, nw)[:flow_direction_selection_ne_p] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        con(wm, nw)[:flow_direction_selection_ne_n] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    con(wm, nw, :flow_direction_selection_ne_p)[a] = Dict{Int, JuMP.ConstraintRef}()
    con(wm, nw, :flow_direction_selection_ne_n)[a] = Dict{Int, JuMP.ConstraintRef}()

    pipe_resistances = ref(wm, nw, :resistance, a)
    constraint_flow_direction_selection_ne(wm, nw, a, pipe_resistances)
end

function constraint_potential_loss_ub_pipe_ne(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    alpha = ref(wm, nw, :alpha)
    len = ref(wm, nw, :pipes, a)["length"]
    pipe_resistances = ref(wm, nw, :resistance, a)

    constraint_potential_loss_ub_pipe_ne(wm, nw, a, alpha, len, pipe_resistances)
end


function constraint_resistance_selection_ne(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    pipe_resistances = ref(wm, nw, :resistance, a)

    constraint_resistance_selection_ne(wm, nw, a, pipe_resistances)
end



### Check Valve Constraints ###
function constraint_check_valve(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, nw), :check_valve_1)
        con(wm, nw)[:check_valve_1] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, nw)[:check_valve_2] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, nw)[:check_valve_3] = Dict{Int, JuMP.ConstraintRef}()
    end

    f_id = ref(wm, nw, :links, a)["f_id"]
    t_id = ref(wm, nw, :links, a)["t_id"]

    constraint_check_valve(wm, nw, a, f_id, t_id)
end

function constraint_potential_loss_check_valve(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, nw), :potential_loss)
        con(wm, nw)[:potential_loss] = Dict{Int, JuMP.ConstraintRef}()
    end

    f_id = ref(wm, nw, :links, a)["f_id"]
    t_id = ref(wm, nw, :links, a)["t_id"]

    len = ref(wm, nw, :pipes, a)["length"]
    r_min = minimum(ref(wm, nw, :resistance, a))

    constraint_potential_loss_check_valve(wm, nw, a, f_id, t_id, len, r_min)
end


### Pump Constraints ###
function constraint_potential_loss_pump(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw, kwargs...)
    f_id = ref(wm, nw, :pumps, a)["f_id"]
    t_id = ref(wm, nw, :pumps, a)["t_id"]

    constraint_potential_loss_pump(wm, nw, a, f_id, t_id)
    constraint_head_gain_pump_quadratic_fit(wm, a; nw=nw, kwargs...)
end

function constraint_head_gain_pump_quadratic_fit(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, nw), :head_gain)
        con(wm, nw)[:head_gain] = Dict{Int, JuMP.ConstraintRef}()
    end

    pump_curve = ref(wm, nw, :pumps, a)["pump_curve"]
    A, B, C = get_function_from_pump_curve(pump_curve)
    constraint_head_gain_pump_quadratic_fit(wm, nw, a, A, B, C)
end
