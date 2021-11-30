#########################################################################
# This file defines commonly-used constraints for water systems models. #
#########################################################################


"""
    constraint_flow_conservation(
        wm, n, i, pipe_fr, pipe_to, des_pipe_fr, des_pipe_to, pump_fr, pump_to,
        regulator_fr, regulator_to, short_pipe_fr, short_pipe_to, valve_fr, valve_to,
        reservoirs, tanks, dispatachable_demands, fixed_demand)

Adds a constraint that ensures flow conservation at a node in the network.
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
    q_pipe, q_des_pipe = var(wm, n, :q_pipe), var(wm, n, :q_des_pipe)
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


"""
    constraint_tank_volume_fixed(wm, n, i, V_0)

Adds a constraint that ensures the volume of a tank at some time step is fixed. Here, `wm`
is the WaterModels object, `n` is the index of a subnetwork within a multinetwork, `i` is
the index of the tank, and `V_0` is the fixed volume of the tank that is desired.
"""
function constraint_tank_volume_fixed(wm::AbstractWaterModel, n::Int, i::Int, V_0::Float64, time_step::Float64, min_vol::Float64, max_vol::Float64)
    V, q = var(wm, n, :V, i), var(wm, n, :q_tank, i)
    c_1 = JuMP.@constraint(wm.model, V == V_0)
    c_2 = JuMP.@constraint(wm.model, min_vol <= V_0 - q * time_step)
    c_3 = JuMP.@constraint(wm.model, V_0 - q * time_step <= max_vol)
    append!(con(wm, n, :tank_volume)[i], [c_1, c_2, c_3])
end

function constraint_tank_volume_fixed_last(wm::AbstractWaterModel, n::Int, i::Int, V_0::Float64)
    V = var(wm, n, :V, i)
    c = JuMP.@constraint(wm.model, V >= V_0)
    con(wm, n, :tank_volume_last)[i] = c
end

function constraint_des_pipe_selection(wm::AbstractWaterModel, n::Int, k::Int, node_fr::Int, node_to::Int, des_pipes::Array{Int64,1})
    z_des_pipe = var(wm, n, :z_des_pipe)
    c = JuMP.@constraint(wm.model, sum(z_des_pipe[a] for a in des_pipes) == 1.0)
    append!(con(wm, n, :des_pipe_selection)[k], [c])
end


"""
    constraint_tank_volume(wm, n_1, n_2, i, time_step)

Adds a constraint that integrates the volume of a tank forward in time. Here, `wm` is the
WaterModels object, `n_1` is the index of a subnetwork within a multinetwork, `n_2` is the
index of another subnetwork forward in time, relative to `n_1`, i is the index of the tank,
and time_step is the time step (in seconds) of the interval from network `n_1` to `n_2`.
"""
function constraint_tank_volume(wm::AbstractWaterModel, n_1::Int, n_2::Int, i::Int, time_step::Float64)
    q_tank = var(wm, n_1, :q_tank, i) # Tank outflow.
    V_1, V_2 = var(wm, n_1, :V, i), var(wm, n_2, :V, i)
    c = JuMP.@constraint(wm.model, V_1 - V_2 == q_tank * time_step)
    append!(con(wm, n_2, :tank_volume)[i], [c])
 end


"""
    constraint_tank_volume_recovery(wm, i, n_1, n_f)

Adds a constraint that ensures the volume of a tank at the end of the time horizon is
greater than or equal to the volume of the tank at the beginning of the time horizon. Here,
`wm` is the WaterModels object, `n_1` is the index of the first subnetwork within a
multinetwork, `n_f` is the index of the final subnetwork, and i is the index of the tank.
"""
function constraint_tank_volume_recovery(wm::AbstractWaterModel, i::Int, n_1::Int, n_f::Int)    
    _initialize_con_dict(wm, :tank_volume_recovery, nw = n_f)
    V_1, V_f = var(wm, n_1, :V, i), var(wm, n_f, :V, i)
    c = JuMP.@constraint(wm.model, V_1 <= V_f)
    con(wm, n_f, :tank_volume_recovery)[i] = c

    # Update the nodal elevation data.
    tank_1, tank_f = ref(wm, n_1, :tank, i), ref(wm, n_f, :tank, i)
    node_1, node_f = ref(wm, n_1, :node, tank_1["node"]), ref(wm, n_f, :node, tank_f["node"])
    node_f["head_min"] = max(node_f["head_min"], node_1["head_min"])
end


function constraint_on_off_pump_power_best_efficiency(
    wm::AbstractWaterModel, n::Int, a::Int,
    density::Float64, gravity::Float64, q_min_forward::Float64)
    # Gather pump flow, power, and status variables.
    q, P, z = var(wm, n, :q_pump, a), var(wm, n, :P_pump, a), var(wm, n, :z_pump, a)

    # Get pump best efficiency data required for construction.
    flow_bep = _calc_pump_best_efficiency_flow(ref(wm, n, :pump, a))
    power_bep = _calc_pump_best_efficiency_power(ref(wm, n, :pump, a), density, gravity)

    # Compute the linear expression used to calculate power.
    P_expr = power_bep * (inv(3.0) * q * inv(flow_bep) + z * 2.0 * inv(3.0))

    # Add constraint equating power with respect to the power curve.
    c = JuMP.@constraint(wm.model, P_expr == P)

    # Append the :on_off_pump_power constraint array.
    append!(con(wm, n, :on_off_pump_power)[a], [c])
end


function constraint_on_off_pump_power_custom(wm::AbstractWaterModel, n::Int, a::Int, power_fixed::Float64, power_variable::Float64)
    # Gather pump flow, power, and status variables.
    q, P, z = var(wm, n, :q_pump, a), var(wm, n, :P_pump, a), var(wm, n, :z_pump, a)
    
    # Add constraint equating power with respect to the linear power curve.
    lhs = power_fixed * z + power_variable * q
    scalar = _get_scaling_factor(vcat(lhs.terms.vals, [1.0]))
    c = JuMP.@constraint(wm.model, scalar * lhs == scalar * P)
 
    # Append the :on_off_pump_power constraint array.
    append!(con(wm, n, :on_off_pump_power)[a], [c])
end


function constraint_on_off_pump_group(wm::AbstractWaterModel, n::Int, k::Int, pump_indices::Set{Int64})
    pump_indices_sorted = sort(collect(pump_indices))

    for i in 1:length(pump_indices_sorted[1:end-1])
        # Add lexicographic constraint for pump statuses.
        z_pump_1 = var(wm, n, :z_pump, pump_indices_sorted[i])
        z_pump_2 = var(wm, n, :z_pump, pump_indices_sorted[i+1])
        c = JuMP.@constraint(wm.model, z_pump_1 >= z_pump_2)
        append!(con(wm, n, :on_off_pump_group)[k], [c])
    end
end


function constraint_on_off_pump_switch(wm::AbstractWaterModel, a::Int, network_ids::Array{Int64, 1}, max_switches::Int64)
    z_switch_on_sum = sum(var(wm, n, :z_switch_on_pump, a) for n in network_ids)
    c = JuMP.@constraint(wm.model, z_switch_on_sum <= max_switches)
    append!(con(wm, network_ids[end], :on_off_pump_switch)[a], [c])
end


function constraint_pump_switch_on(wm::AbstractWaterModel, a::Int, n_1::Int, n_2::Int, nws_active::Array{Int64, 1})
    z_1, z_2 = var(wm, n_1, :z_pump, a), var(wm, n_2, :z_pump, a)
    z_switch_on = var(wm, n_2, :z_switch_on_pump, a)
    c_1 = JuMP.@constraint(wm.model, z_switch_on >= z_2 - z_1)
    append!(con(wm, n_2, :pump_switch_on)[a], [c_1])

    for nw_active in nws_active
        z_nw = var(wm, nw_active, :z_pump, a)
        c_2 = JuMP.@constraint(wm.model, z_switch_on <= z_nw)
        append!(con(wm, n_2, :pump_switch_on)[a], [c_2])
    end
end


function constraint_pump_switch_off(wm::AbstractWaterModel, a::Int, n_1::Int, n_2::Int, nws_inactive::Array{Int64, 1})
    z_1, z_2 = var(wm, n_1, :z_pump, a), var(wm, n_2, :z_pump, a)
    z_switch_off = var(wm, n_2, :z_switch_off_pump, a)
    c_1 = JuMP.@constraint(wm.model, z_switch_off >= z_1 - z_2)
    append!(con(wm, n_2, :pump_switch_off)[a], [c_1])

    for nw_inactive in nws_inactive
        z_nw = var(wm, nw_inactive, :z_pump, a)
        c_2 = JuMP.@constraint(wm.model, z_nw <= 1.0 - z_switch_off)
        append!(con(wm, n_2, :pump_switch_off)[a], [c_2])
    end
end


"Try to determine a scaling factor that centers values around one."
function _get_scaling_factor(values::Vector{Float64})::Float64
    values_filtered = values[values .!= 0.0]
    mean_log10_value = Statistics.mean(log10.(abs.(values_filtered)))
    return 10^(-mean_log10_value)
end
