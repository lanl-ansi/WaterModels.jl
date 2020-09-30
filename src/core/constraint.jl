#########################################################################
# This file defines commonly-used constraints for water systems models. #
#########################################################################


"""
    constraint_reservoir_head(wm, n, i, head)
"""
function constraint_reservoir_head(wm::AbstractWaterModel, n::Int, i::Int, head::Float64)
    h = var(wm, n, :h, i)
    c = JuMP.@constraint(wm.model, h == head)
    con(wm, n, :reservoir_head)[i] = c
end


"""
    constraint_flow_conservation(
        wm, n, i, pipe_fr, pipe_to, des_pipe_fr, des_pipe_to, pump_fr, pump_to,
        regulator_fr, regulator_to, short_pipe_fr, short_pipe_to, valve_fr, valve_to,
        reservoirs, tanks, dispatachable_demands, fixed_demand)
"""
function constraint_flow_conservation(
    wm::AbstractWaterModel, n::Int, i::Int, pipe_fr::Array{Int64,1},
    pipe_to::Array{Int64,1}, des_pipe_fr::Array{Int64,1}, des_pipe_to::Array{Int64,1},
    pump_fr::Array{Int64,1}, pump_to::Array{Int64,1}, regulator_fr::Array{Int64,1},
    regulator_to::Array{Int64,1}, short_pipe_fr::Array{Int64,1},
    short_pipe_to::Array{Int64,1}, valve_fr::Array{Int64,1}, valve_to::Array{Int64,1},
    reservoirs::Array{Int64,1}, tanks::Array{Int64,1}, dispatchable_demands::Array{Int64,1},
    fixed_demand::Float64)
    # Collect flow variable references per component.
    q_pipe, q_des_pipe = var(wm, n, :q_pipe), var(wm, n, :q_des_pipe_sum)
    q_pump, q_regulator = var(wm, n, :q_pump), var(wm, n, :q_regulator)
    q_short_pipe, q_valve = var(wm, n, :q_short_pipe), var(wm, n, :q_valve)
    q_reservoir, q_tank = var(wm, n, :q_reservoir), var(wm, n, :q_tank)
    q_demand = var(wm, n, :q_demand)

    # Add the flow conservation constraint.
    con(wm, n, :flow_conservation)[i] = JuMP.@constraint(wm.model, -
         sum(q_pipe[a] for a in pipe_fr) + sum(q_pipe[a] for a in pipe_to) -
         sum(q_des_pipe[a] for a in des_pipe_fr) + sum(q_des_pipe[a] for a in des_pipe_to) -
         sum(q_pump[a] for a in pump_fr) + sum(q_pump[a] for a in pump_to) -
         sum(q_regulator[a] for a in regulator_fr) +
         sum(q_regulator[a] for a in regulator_to) -
         sum(q_short_pipe[a] for a in short_pipe_fr) +
         sum(q_short_pipe[a] for a in short_pipe_to) -
         sum(q_valve[a] for a in valve_fr) + sum(q_valve[a] for a in valve_to) == -
         sum(q_reservoir[k] for k in reservoirs) - sum(q_tank[k] for k in tanks) +
         sum(q_demand[k] for k in dispatchable_demands) + fixed_demand)
end


function constraint_pump_flow(wm::AbstractWaterModel, n::Int, a::Int, q_min_active::Float64)
    # Get pump status variable.
    q, z = var(wm, n, :q_pump), var(wm, n, :z_pump, a)

    # If the pump is inactive, flow must be zero.
    q_ub = JuMP.upper_bound(q)
    c_1 = JuMP.@constraint(wm.model, q >= q_min_active * z)
    c_2 = JuMP.@constraint(wm.model, q <= q_ub * z)

    # Append the constraint array.
    append!(con(wm, n, :on_off_pump_flow, a), [c_1, c_2])
end


function constraint_pump_head(wm::AbstractWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int)
    # Get head variables for from and to nodes.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # Get pump status variable.
    g, z = var(wm, n, :g_pump), var(wm, n, :z_pump, a)

    # If the pump is off, decouple the head difference relationship.
    dh_ub = JuMP.upper_bound(h_i) - JuMP.lower_bound(h_j)
    dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    c_1 = JuMP.@constraint(wm.model, h_i - h_j <= dh_ub * (1.0 - z))
    c_2 = JuMP.@constraint(wm.model, h_i - h_j >= g + dh_lb * (1.0 - z))
    c_3 = JuMP.@constraint(wm.model, h_i - h_j <= g)

    # Append the constraint array.
    append!(con(wm, n, :pump_head, a), [c_1, c_2, c_3])
end


function constraint_tank_volume_initial(wm::AbstractWaterModel, n::Int, i::Int, V_0::Float64)
    V = var(wm, n, :V, i)
    c = JuMP.@constraint(wm.model, V == V_0)
    con(wm, n, :tank_volume)[i] = c
end


"""
    constraint_tank_volume(wm, n_1, n_2, i, time_step)

Adds a constraint that integrates the volume of a tank forward in time. Here, `wm` is the
WaterModels object, `n_1` is the index of a subnetwork within a multinetwork, `n_2` is the
index of another subnetwork forward in time, relative to `n_1`, i is the index of the tank,
and time_step is the time step (in seconds) between networks `n_1` and `n_2`.
"""
function constraint_tank_volume(wm::AbstractWaterModel, n_1::Int, n_2::Int, i::Int, time_step::Float64)
    qt = var(wm, n_1, :q_tank, i) # Tank outflow.
    V_1, V_2 = var(wm, n_1, :V, i), var(wm, n_2, :V, i)
    c = JuMP.@constraint(wm.model, V_1 - time_step * qt == V_2)
    con(wm, n_2, :tank_volume)[i] = c
end


function constraint_tank_volume_recovery(wm::AbstractWaterModel, i::Int, n_1::Int, n_f::Int)
    if !ref(wm, nw_f, :tank, i)["dispatchable"]
        _initialize_con_dict(wm, :tank_volume_recovery, nw=n_f)
        V_1, V_f = var(wm, n_1, :V, i), var(wm, n_f, :V, i)
        c = JuMP.@constraint(wm.model, V_f >= V_1)
        con(wm, n_f, :tank_volume_recovery)[i] = c
    end
end
