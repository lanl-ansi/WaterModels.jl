#########################################################################
# This file defines commonly-used constraints for water systems models. #
#########################################################################


function constraint_source_head(wm::AbstractWaterModel, n::Int, i::Int, h_s::Float64)
    c = JuMP.@constraint(wm.model, var(wm, n, :h, i) == h_s)
    con(wm, n, :source_head)[i] = c
end


"""
    constraint_flow_conservation(
        wm, n, i, check_valve_fr, check_valve_to, pipe_fr, pipe_to, pump_fr, pump_to,
        pressure_reducing_valve_fr, pressure_reducing_valve_to, shutoff_valve_fr,
        shutoff_valve_to, reservoirs, tanks, demand)
"""
function constraint_flow_conservation(
    wm::AbstractWaterModel, n::Int, i::Int, check_valve_fr::Array{Int64,1},
    check_valve_to::Array{Int64,1}, pipe_fr::Array{Int64,1}, pipe_to::Array{Int64,1},
    pump_fr::Array{Int64,1}, pump_to::Array{Int64,1},
    pressure_reducing_valve_fr::Array{Int64,1}, pressure_reducing_valve_to::Array{Int64,1},
    shutoff_valve_fr::Array{Int64,1}, shutoff_valve_to::Array{Int64,1},
    reservoirs::Array{Int64,1}, tanks::Array{Int64,1}, demand::Float64)
    # Collect flow variable references per component.
    q_check_valve = var(wm, n, :q_check_valve)
    q_pipe, q_pump = var(wm, n, :q_pipe), var(wm, n, :q_pump)
    q_pressure_reducing_valve = var(wm, n, :q_pressure_reducing_valve)
    q_shutoff_valve = var(wm, n, :q_shutoff_valve)

    # Add the flow conservation constraint.
    con(wm, n, :flow_conservation)[i] = JuMP.@constraint(wm.model,
        - sum(q_check_valve[a] for a in check_valve_fr)
        + sum(q_check_valve[a] for a in check_valve_to)
        - sum(q_pipe[a] for a in pipe_fr) + sum(q_pipe[a] for a in pipe_to)
        - sum(q_pump[a] for a in pump_fr) + sum(q_pump[a] for a in pump_to)
        - sum(q_pressure_reducing_valve[a] for a in pressure_reducing_valve_fr)
        + sum(q_pressure_reducing_valve[a] for a in pressure_reducing_valve_to)
        - sum(q_shutoff_valve[a] for a in shutoff_valve_fr)
        + sum(q_shutoff_valve[a] for a in shutoff_valve_to) ==
        - sum(var(wm, n, :qr, id) for id in reservoirs)
        - sum(var(wm, n, :qt, id) for id in tanks) + demand)
end


function constraint_volume(wm::AbstractWaterModel, n::Int, i::Int, elevation::Float64, surface_area::Float64)
    h, V = var(wm, n, :h, i), var(wm, n, :V, i)
    c = JuMP.@constraint(wm.model, h - elevation == V * inv(surface_area))
    con(wm, n, :volume)[i] = c
end


function constraint_pump_control_initial(wm::AbstractWaterModel, n::Int, a::Int, status::Bool)
    z = var(wm, n, :z_pump, a)
    c = JuMP.@constraint(wm.model, z == status)
end


function constraint_tank_state_initial(wm::AbstractWaterModel, n::Int, i::Int, V_0::Float64, time_step::Float64)
    V = var(wm, n, :V, i)
    c = JuMP.@constraint(wm.model, V == V_0)
    con(wm, n, :tank_state)[i] = c
end


"""
    constraint_tank_state(wm, n_1, n_2, i, time_step)

Adds a constraint that integrates the volume of a tank forward in time. Here, `wm` is the
WaterModels object, `n_1` is the index of a subnetwork within a multinetwork, `n_2` is the
index of another subnetwork forward in time, relative to `n_1`, i is the index of the tank,
and time_step is the time step (in seconds) between networks `n_1` and `n_2`.
"""
function constraint_tank_state(wm::AbstractWaterModel, n_1::Int, n_2::Int, i::Int, time_step::Float64)
    qt = var(wm, n_1, :qt, i) # Tank outflow.
    V_1, V_2 = var(wm, n_1, :V, i), var(wm, n_2, :V, i)
    c = JuMP.@constraint(wm.model, V_1 - time_step * qt == V_2)
    con(wm, n_2, :tank_state)[i] = c
end


function constraint_recover_volume(wm::AbstractWaterModel, i::Int, n_1::Int, n_f::Int)
    _initialize_con_dict(wm, :recover_volume, nw=n_f)
    V_1, V_f = var(wm, n_1, :V, i), var(wm, n_f, :V, i)
    c = JuMP.@constraint(wm.model, V_f >= V_1)
    con(wm, n_f, :recover_volume)[i] = c
end
