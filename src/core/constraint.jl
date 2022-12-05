#########################################################################
# This file defines commonly-used constraints for water systems models. #
#########################################################################


"""
    constraint_flow_conservation(
        wm::AbstractWaterModel,
        n::Int,
        i::Int,
        pipe_fr::Array{Int,1},
        pipe_to::Array{Int,1},
        des_pipe_fr::Array{Int,1},
        des_pipe_to::Array{Int,1},
        pump_fr::Array{Int,1},
        pump_to::Array{Int,1},
        regulator_fr::Array{Int,1},
        regulator_to::Array{Int,1},
        short_pipe_fr::Array{Int,1},
        short_pipe_to::Array{Int,1},
        valve_fr::Array{Int,1},
        valve_to::Array{Int,1},
        reservoirs::Array{Int,1},
        tanks::Array{Int,1},
        dispatchable_demands::Array{Int,1},
        fixed_demand::Float64
    )

Adds a constraint that ensures volumetric flow rate (and thus mass) is
conserved at node `i` and subnetwork (or time) index `n` in the network. Here,
`pipe_fr`, `pipe_to`, etc., are components that are directed _from_ or _to_
node `i`, respectively. Additionally, `reservoirs`, `tanks`, and
`dispatchable_demands` are node-attached components that variably contribute to
nodal inflow and outflow at node `i`. Finally, `fixed_demand` is the fixed
amount of flow demanded at node `i`, which may be negative _or_ nonnegative.
"""
function constraint_flow_conservation(
    wm::AbstractWaterModel, n::Int, i::Int, pipe_fr::Array{Int,1},
    pipe_to::Array{Int,1}, des_pipe_fr::Array{Int,1}, des_pipe_to::Array{Int,1},
    pump_fr::Array{Int,1}, pump_to::Array{Int,1}, regulator_fr::Array{Int,1},
    regulator_to::Array{Int,1}, short_pipe_fr::Array{Int,1},
    short_pipe_to::Array{Int,1}, valve_fr::Array{Int,1}, valve_to::Array{Int,1},
    reservoirs::Array{Int,1}, tanks::Array{Int,1}, dispatchable_demands::Array{Int,1},
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
    constraint_tank_volume_fixed(
        wm::AbstractWaterModel,
        n::Int,
        i::Int,
        V_0::Float64,
        time_step::Float64,
        min_vol::Float64,
        max_vol::Float64
    )

Adds a constraint that ensures the volume of a tank at some time index is fixed.
Here, `wm` is the WaterModels object, `n` is the index of a subnetwork within a
multinetwork, `i` is the index of the tank, and `V_0` is the fixed volume of
the tank that is desired. Also adds constraints that ensure, after integration
of volume forward in time, the new volume will reside between predefined lower
and upper bounds. To that end, `time_step` is the time step between the current
time index and the next, `min_vol` is the minimum volume of water that must be
present in the tank, and `max_vol` is the maximum volume of water in the tank.
"""
function constraint_tank_volume_fixed(wm::AbstractWaterModel, n::Int, i::Int, V_0::Float64, time_step::Float64, min_vol::Float64, max_vol::Float64)
    V, q = var(wm, n, :V, i), var(wm, n, :q_tank, i)
    c_1 = JuMP.@constraint(wm.model, V == V_0)
    c_2 = JuMP.@constraint(wm.model, min_vol <= V_0 - q * time_step)
    c_3 = JuMP.@constraint(wm.model, V_0 - q * time_step <= max_vol)
    append!(con(wm, n, :tank_volume)[i], [c_1, c_2, c_3])
end


"""
    constraint_des_pipe_selection(
        wm::AbstractWaterModel,
        n::Int,
        k::Int,
        node_fr::Int,
        node_to::Int,
        des_pipes::Array{Int,1}
    )

Adds a constraint that ensures exactly one design pipe will be selected per arc
that comprises one or more design pipe possibilities. Here, `wm` is the
WaterModels object, `n` is the index of a subnetwork within a multinetwork, `k`
is the index of the arc comprising multiple design pipes, `node_fr` and
`node_to` are node indices that connect the arc, and `des_pipes` is the vector
of design pipes that reside along the same arc `k`.
"""
function constraint_des_pipe_selection(wm::AbstractWaterModel, n::Int, k::Int, node_fr::Int, node_to::Int, des_pipes::Array{Int,1})
    z_des_pipe = var(wm, n, :z_des_pipe)
    c = JuMP.@constraint(wm.model, sum(z_des_pipe[a] for a in des_pipes) == 1.0)
    append!(con(wm, n, :des_pipe_selection)[k], [c])
end


"""
    constraint_tank_volume(
        wm::AbstractWaterModel,
        n_1::Int,
        n_2::Int,
        i::Int,
        time_step::Float64
    )

Adds a constraint that integrates the volume of a tank forward in time. Here,
`wm` is the WaterModels object, `n_1` is the index of a subnetwork within a
multinetwork, `n_2` is the index of another subnetwork forward in time,
relative to `n_1`, i is the index of the tank, and `time_step` is the time step
of the interval from network `n_1` to `n_2`.
"""
function constraint_tank_volume(wm::AbstractWaterModel, n_1::Int, n_2::Int, i::Int, time_step::Float64)
    q_tank = var(wm, n_1, :q_tank, i) # Tank _outflow_.
    V_1, V_2 = var(wm, n_1, :V, i), var(wm, n_2, :V, i)
    c = JuMP.@constraint(wm.model, V_1 - V_2 == q_tank * time_step)
    append!(con(wm, n_2, :tank_volume)[i], [c])
 end


"""
    constraint_tank_volume_recovery(
        wm::AbstractWaterModel,
        i::Int,
        n_1::Int,
        n_f::Int
    )

Adds a constraint that ensures the volume of a tank at the end of the time
horizon is greater than or equal to the volume of the tank at the beginning
of the time horizon. Here, `wm` is the WaterModels object, `i` is the index of
the tank, `n_1` is the index of the first subnetwork within a multinetwork, and
`n_f` is the index of the final subnetwork. Also updates the minimum head in
the `ref` dictionary for node `i` at the final subnetwork index `n_f` to take
the minimum head value at `n_1` if it is larger, as implied by this constraint.
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


"""
    constraint_on_off_pump_power_best_efficiency(
        wm::AbstractWaterModel,
        n::Int,
        a::Int,
        density::Float64,
        gravity::Float64,
        q_min_forward::Float64
    )

Adds a constraint that expresses pump power as a linear expression of flow,
where coefficients of this linear expression are computed using an assumption
that the pump will be operating at its best efficiency point. Here, `wm` is the
WaterModels object, `n` is the index of the subnetwork within a multinetwork,
`a` is the index of the pump, `density` is the density of water, `gravity` is
acceleration due to gravity, and `q_min_forward` is the minimum amount of flow
transported through the pump when it is active (and flow is directed forward).
"""
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


"""
    constraint_on_off_pump_power_custom(
        wm::AbstractWaterModel,
        n::Int,
        a::Int,
        power_fixed::Float64,
        power_variable::Float64
    )

Adds a constraint that expresses pump power as a linear expression of flow,
where coefficients of this linear expression are given by `power_fixed` and
`power_variable`. Here, `wm` is the WaterModels object, `n` is the index of the
subnetwork within a multinetwork, `a` is the index of the pump, `power_fixed`
is the amount of power consumed by the pump when it is active and transporting
no flow (i.e., the point of intersection at zero flow on the vertical axis of a
power-versus-flow curve), and `power_variable` is the amount of power consumed
by the pump per unit flow.
"""
function constraint_on_off_pump_power_custom(wm::AbstractWaterModel, n::Int, a::Int, power_fixed::Float64, power_variable::Float64)
    # Gather pump flow, power, and status variables.
    q, P, z = var(wm, n, :q_pump, a), var(wm, n, :P_pump, a), var(wm, n, :z_pump, a)
    
    # Add constraint equating power with respect to the linear power curve.
    c = JuMP.@constraint(wm.model, power_fixed * z + power_variable * q == P)

    # Append the :on_off_pump_power constraint array.
    append!(con(wm, n, :on_off_pump_power)[a], [c])
end


"""
    constraint_on_off_pump_group(
        wm::AbstractWaterModel,
        n::Int,
        k::Int,
        pump_indices::Set{Int}
    )

Adds a constraint that imposes symmetry-breaking lexicographic sorting of pump
activation statuses on groups of identical pumps operating in parallel along
the same arc of the network. Reduces the combinatorial complexity of problems
involving multiple pumps operating in parallel. Here, `wm` is the WaterModels
object, `n` is the index of the subnetwork within a multinetwork, `k` is the
index of the arc that includes the group of identical pumps, and `pump_indices`
is the set of pump indices for pumps operating in parallel along arc `k`.
"""
function constraint_on_off_pump_group(wm::AbstractWaterModel, n::Int, k::Int, pump_indices::Set{Int})
    pump_indices_sorted = sort(collect(pump_indices))

    for i in 1:length(pump_indices_sorted[1:end-1])
        # Add lexicographic constraint for pump statuses.
        z_pump_1 = var(wm, n, :z_pump, pump_indices_sorted[i])
        z_pump_2 = var(wm, n, :z_pump, pump_indices_sorted[i+1])
        c = JuMP.@constraint(wm.model, z_pump_1 >= z_pump_2)
        append!(con(wm, n, :on_off_pump_group)[k], [c])
    end
end


"""
    constraint_on_off_pump_switch(
        wm::AbstractWaterModel,
        a::Int,
        network_ids::Array{Int,1},
        max_switches::Int
    )

Adds a constraint that limits the number of times a pump can be switched from
off to on in a multiperiod pump scheduling problem. Here, `wm` is the
WaterModels object, `a` is the index of the pump, `network_ids` are the network
(time) indices used in the summation that limits the number of switches, and
`max_switches` is the number of maximum switches permitted for the pump over
the set of network indices.
"""
function constraint_on_off_pump_switch(wm::AbstractWaterModel, a::Int, network_ids::Array{Int, 1}, max_switches::Int)
    z_switch_on_sum = sum(var(wm, n, :z_switch_on_pump, a) for n in network_ids)
    c = JuMP.@constraint(wm.model, z_switch_on_sum <= max_switches)
    append!(con(wm, network_ids[end], :on_off_pump_switch)[a], [c])
end


"""
    constraint_pump_switch_on(
        wm::AbstractWaterModel,
        a::Int,
        n_1::Int,
        n_2::Int,
        nws_active::Array{Int,1}
    )

Adds a constraint that models the switching of a pump from _off_ to _on_ and
constrains its operation, if switched on, to remain on for at least some amount
of time. Here, `wm` is the WaterModels object; `a` is the index of the pump;
`n_1` is the first time index modeled by the constraint; `n_2` is the adjacent
(next) time index modeled by the constraint, which permits limiting the change
in pump status (i.e., from off to on); and `nws_active` are the subnetwork
(time) indices where the pump must be constrained to on _if_ the pump is indeed
switched from off to on between time indices `n_1` and `n_2`.
"""
function constraint_pump_switch_on(wm::AbstractWaterModel, a::Int, n_1::Int, n_2::Int, nws_active::Array{Int, 1})
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


"""
    constraint_pump_switch_off(
        wm::AbstractWaterModel,
        a::Int,
        n_1::Int,
        n_2::Int,
        nws_inactive::Array{Int,1}
    )

Adds a constraint that models the switching of a pump from _on_ to _off_ and
constrains non-operation, if switched off, to remain off for at least some
amount of time. Here, `wm` is the WaterModels object; `a` is the index of the
pump; `n_1` is the first time index modeled by the constraint; `n_2` is the
adjacent (next) time index modeled by the constraint, which permits limiting
the change in pump status (i.e., from on to off); and `nws_inactive` are the
subnetwork (time) indices where the pump must be constrained to off _if_ the
pump is indeed switched from on to off between time indices `n_1` and `n_2`.
"""
 function constraint_pump_switch_off(wm::AbstractWaterModel, a::Int, n_1::Int, n_2::Int, nws_inactive::Array{Int, 1})
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
