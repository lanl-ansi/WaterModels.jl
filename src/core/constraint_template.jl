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

    # TBD
    # for tid in ref(wm, n, :node_tanks, i)
    #     tank = ref(wm, n, :tanks, tid)
    #     # TODO add tank vars as loads
    # end

    node_demands = Dict(k => ref(wm, nw, :junctions, k, "demand") for k in node_junctions)

    constraint_flow_conservation(wm, nw, i, node_arcs_fr, node_arcs_to, node_reservoirs, node_demands)
end


function constraint_sink_flow(wm::GenericWaterModel, i::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, n), :directed_sink_flow)
        con(wm, n)[:directed_sink_flow] = Dict{Int, JuMP.ConstraintRef}()
    end

    constraint_sink_flow(wm, nw, i)
end

function constraint_source_flow(wm::GenericWaterModel, i::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, n), :directed_source_flow)
        con(wm, n)[:directed_source_flow] = Dict{Int, JuMP.ConstraintRef}()
    end

    constraint_sink_flow(wm, nw, i)
end



### Junction Constraints ###



### Reservoir Constraints ###



### Tank Constraints ###
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

    if haskey(ref(wm, nw), :time_series)
        time_step = ref(wm, nw, :time_series)["time_step"]
    else
        Memento.error(_LOGGER, "Tank states cannot be controlled outside of a time series.")
    end

    constraint_tank_state(wm, nw_1, nw_2, i, time_step)
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
    constraint_potential_loss_pipe(wm, nw, a)

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
    constraint_potential_loss_pipe_ne(wm, nw, a)

    constraint_head_difference(wm, a; nw=nw, kwargs...)
    constraint_flow_direction_selection_ne(wm, a; nw=nw, kwargs...)
    constraint_potential_loss_ub_pipe_ne(wm, a; nw=nw, kwargs...)
end

function constraint_flow_direction_selection_ne(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    constraint_flow_direction_selection_ne(wm, nw, a)
end

function constraint_potential_loss_ub_pipe_ne(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    constraint_potential_loss_ub_pipe_ne(wm, nw, a)
end


function constraint_resistance_selection_ne(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    constraint_resistance_selection_ne(wm, nw, a)
end



### Check Valve Constraints ###
function constraint_check_valve(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, nw), :check_valve_1)
        con(wm, nw)[:check_valve_1] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, nw)[:check_valve_2] = Dict{Int, JuMP.ConstraintRef}()
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


