########################################################################
# This file defines commonly-used variables for water systems models.
########################################################################


"Sets the start value for a given variable."
function comp_start_value(comp::Dict{String,<:Any}, key::String, default::Float64=0.0)
    return get(comp, key, default)
end


"Sets the start value for a given variable."
function comp_start_value(comp::Dict{String,<:Any}, key_1::String, key_2::Int, default=0.0)
    return key_1 in keys(comp) ? get(get(comp, key_1, default), key_2, default) : default
end


### Variables related to nodal components. ###

"""
    variable_head(
        wm::AbstractWaterModel;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true
    )

Instantiates bounded (by default) or unbounded total hydraulic head (or head)
variables for nodes in the network at subnetwork (or time) index `nw`, i.e.,
`h[i]` for `i` in `node`. Also instantiates JuMP expressions that are affine
expressions of each head, i.e., pressure, which is computed as `p[i] = h[i] -
elevation[i]`, and tank volume (at each tank attached to a node), which is
computed as `V[k] = 0.25 * pi * diameter[k]^2 * (h[i] - elevation[i])`, where
`k` is the index of the tank; `diameter[k]` is the cross-sectional diameter of
the (cylindrical) tank; and `i` is the index of the node to which the tank is
attached.
"""
function variable_head(wm::AbstractWaterModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # Initialize variables for total hydraulic head.
    h = var(wm, nw)[:h] = JuMP.@variable(wm.model,
        [i in ids(wm, nw, :node)], base_name = "$(nw)_h",
        start = comp_start_value(ref(wm, nw, :node, i), "h_start"))

    # Initialize variables for head differences across design pipe arcs.
    dh_des_pipe = var(wm, nw)[:dh_des_pipe] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :des_pipe)], base_name = "$(nw)_dh",
        start = comp_start_value(ref(wm, nw, :des_pipe, a), "dh_start"))

    if bounded
        for (i, node) in ref(wm, nw, :node)
            # Set the lower and upper bounds for each head.
            JuMP.set_lower_bound(h[i], node["head_min"])
            JuMP.set_upper_bound(h[i], node["head_max"])

            # Set start value for the head variable with possibly better data.
            h_mid = 0.5 * (node["head_max"] + node["head_min"])
            h_start = comp_start_value(node, "h_start", h_mid)
            JuMP.set_start_value(h[i], h_start)
        end

        for (a, des_pipe) in ref(wm, nw, :des_pipe)
            # Get the nodes that are connected by the design pipe.
            node_fr = ref(wm, nw, :node, des_pipe["node_fr"])
            node_to = ref(wm, nw, :node, des_pipe["node_to"])

            # Set the lower and upper bounds for the design pipe head difference.
            JuMP.set_lower_bound(dh_des_pipe[a], node_fr["head_min"] - node_to["head_max"])
            JuMP.set_upper_bound(dh_des_pipe[a], node_fr["head_max"] - node_to["head_min"])
        end
    end

    # Initialize an entry to the solution component dictionary for heads.
    report && sol_component_value(wm, nw, :node, :h, ids(wm, nw, :node), h)

    # Create expressions that calculate pressures from head.
    p = var(wm, nw)[:p] = JuMP.@expression(wm.model,
        [i in ids(wm, nw, :node)], var(wm, nw, :h, i) -
        ref(wm, nw, :node, i)["elevation"])

    # Initialize an entry to the solution component dictionary for pressures.
    report && sol_component_value(wm, nw, :node, :p, ids(wm, nw, :node), p)

    # Create an expression that maps total hydraulic head to volume for tanks.
    V = var(wm, nw)[:V] = JuMP.@expression(wm.model, [i in ids(wm, nw, :tank)],
        (var(wm, nw, :h, ref(wm, nw, :tank, i)["node"]) -
        ref(wm, nw, :node, ref(wm, nw, :tank, i)["node"])["elevation"]) *
        0.25 * pi * ref(wm, nw, :tank, i)["diameter"]^2)

    # Initialize an entry to the solution component dictionary for volumes.
    report && sol_component_value(wm, nw, :tank, :V, ids(wm, nw, :tank), V)
end


"""
    variable_pump_head_gain(
        wm::AbstractWaterModel;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true
    )

Creates bounded (by default) or unbounded head gain variables for pumps in
in the network at subnetwork (or time) index `nw`, i.e., `g_pump[a]` for `a`
in `pump`. Note that these variables are always nonnegative, since for each
pump, there will never be negative head gain (i.e., pumps only increase head).
Head gain is always directed from "node_fr" (i.e., the tail node of the arc) to
"node_to" (i.e., the head node of the arc).
"""
function variable_pump_head_gain(wm::AbstractWaterModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # Initialize variables for total hydraulic head gain from a pump.
    g = var(wm, nw)[:g_pump] = JuMP.@variable(wm.model, [a in ids(wm, nw, :pump)],
        base_name="$(nw)_g_pump", lower_bound=0.0, # Pump gain is nonnegative.
        start=comp_start_value(ref(wm, nw, :pump, a), "g_pump_start", 1.0e-6))

    if bounded
        for (a, pump) in ref(wm, nw, :pump)
            # Get the nodes that are connected by the pump.
            node_fr = ref(wm, nw, :node, pump["node_fr"])
            node_to = ref(wm, nw, :node, pump["node_to"])

            # Set the upper bound for the head gain variable.
            head_gain_max = _calc_pump_head_gain_max(pump, node_fr, node_to)
            JuMP.set_upper_bound(g[a], head_gain_max)

            # Set the start value for the head gain variable with possibly better data.
            g_mid = 0.5 * head_gain_max # Midpoint of possible head gains.
            g_start = comp_start_value(pump, "g_pump_start", g_mid)
            JuMP.set_start_value(g[a], g_start)
        end
    end

    # Initialize an entry to the solution component dictionary for head gains.
    report && sol_component_value(wm, nw, :pump, :g, ids(wm, nw, :pump), g)
end


"""
    variable_pump_power(
        wm::AbstractWaterModel;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true
    )

Creates bounded (by default) or unbounded power consumption variables for pumps
in the network at subnetwork (or time) index `nw`, i.e., `P[a]` for `a` in
`pump`. Note that these variables are always nonnegative since each pump only
_consumes_ power. Additionally, two sets of JuMP expressions are also derived
from these power variables: energy consumption, i.e., `E[a]` for `a` in `pump`
and cost of pump operation across a time step, i.e., `c[a]` for `a` in `pump`.
"""
function variable_pump_power(wm::AbstractWaterModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # Compute scaled density and gravity.
    wm_data = get_wm_data(wm.data)
    base_mass = wm_data["per_unit"] ? wm_data["base_mass"] : 1.0
    base_length = wm_data["per_unit"] ? wm_data["base_length"] : 1.0
    base_time = wm_data["per_unit"] ? wm_data["base_time"] : 1.0

    rho_s = _calc_scaled_density(base_mass, base_length)
    g_s = _calc_scaled_gravity(base_length, base_time)

    # Initialize variables for the power utilization of a pump.
    P = var(wm, nw)[:P_pump] = JuMP.@variable(wm.model, [a in ids(wm, nw, :pump)],
        base_name = "$(nw)_P_pump", lower_bound = 0.0, # Pump power is nonnegative.
        start = comp_start_value(ref(wm, nw, :pump, a), "P_pump_start", 0.0))

    if bounded
        for (a, pump) in ref(wm, nw, :pump)
            # Get the nodes that are connected by the pump.
            node_fr = ref(wm, nw, :node, pump["node_fr"])
            node_to = ref(wm, nw, :node, pump["node_to"])
    
            # Set the upper bound for the power variable.
            P_max = _calc_pump_power_max(pump, node_fr, node_to, rho_s, g_s)
            JuMP.set_upper_bound(P[a], P_max)

            # Set the start value for the scaled power variable.
            P_start = comp_start_value(pump, "P_pump_start", 0.5 * P_max)
            JuMP.set_start_value(P[a], P_start)
        end
    end

    # Create expressions to compute the unscaled pump energies.
    E = var(wm, nw)[:E_pump] = JuMP.@expression(
        wm.model, [a in ids(wm, nw, :pump)],
        P[a] * ref(wm, nw, :time_step))

    # Create expressions to compute the unscaled pump costs.
    c = var(wm, nw)[:c_pump] = JuMP.@expression(
        wm.model, [a in ids(wm, nw, :pump)], E[a] *
        ref(wm, nw, :pump, a)["energy_price"])

    if report
        # Initialize entries to the solution component dictionary for expressions.
        sol_component_value(wm, nw, :pump, :P, ids(wm, nw, :pump), P)
        sol_component_value(wm, nw, :pump, :E, ids(wm, nw, :pump), E)
        sol_component_value(wm, nw, :pump, :c, ids(wm, nw, :pump), c)
    end
end


"""
    variable_reservoir_flow(
        wm::AbstractWaterModel;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true
    )

Creates bounded (by default) or unbounded outgoing volumetric flow rate
variables for reservoirs in the network at subnetwork (or time) index `nw`,
i.e., `q_reservoir[i]` for `i` in `reservoir`. Note that these variables are
always nonnegative, since for each reservoir, there will never be incoming
flow.
"""
function variable_reservoir_flow(wm::AbstractWaterModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    q_reservoir = var(wm, nw)[:q_reservoir] = JuMP.@variable(wm.model,
        [i in ids(wm, nw, :reservoir)], lower_bound=0.0, base_name="$(nw)_q_reservoir",
        start=comp_start_value(ref(wm, nw, :reservoir, i), "q_reservoir_start"))

    if bounded
        for (i, reservoir) in ref(wm, nw, :reservoir)
            flow_min, flow_max = 0.0, 0.0

            for name in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
                edges_fr = ref(wm, nw, Symbol(name * "_fr"), reservoir["node"])
                edges_to = ref(wm, nw, Symbol(name * "_to"), reservoir["node"])

                if length(edges_fr) > 0
                    flow_min += sum(ref(wm, nw, Symbol(name), a)["flow_min"] for a in edges_fr)
                    flow_max += sum(ref(wm, nw, Symbol(name), a)["flow_max"] for a in edges_fr)
                end

                if length(edges_to) > 0
                    flow_min -= sum(ref(wm, nw, Symbol(name), a)["flow_max"] for a in edges_to)
                    flow_max -= sum(ref(wm, nw, Symbol(name), a)["flow_min"] for a in edges_to)
                end
            end

            flow_min = max(0.0, flow_min)
            JuMP.set_lower_bound(q_reservoir[i], flow_min)
            JuMP.set_upper_bound(q_reservoir[i], flow_max)

            flow_mid = 0.5 * (flow_max + flow_min)
            q_start = comp_start_value(reservoir, "q_reservoir_start", flow_mid)
            JuMP.set_start_value(q_reservoir[i], q_start)
        end
    end

    report && sol_component_value(wm, nw, :reservoir, :q, ids(wm, nw, :reservoir), q_reservoir)
end


"""
    variable_demand_flow(
        wm::AbstractWaterModel;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true
    )

Creates unbounded volumetric flow rate demand variables for all _dispatchable_
demands in the network at subnetwork (or time) index `nw`, i.e., `q_demand[i]`
for `i` in `dispatchable_demand`.
"""
function variable_demand_flow(wm::AbstractWaterModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    q_demand = var(wm, nw)[:q_demand] = JuMP.@variable(wm.model,
        [i in ids(wm, nw, :dispatchable_demand)], base_name="$(nw)_q_demand",
        start=comp_start_value(ref(wm, nw, :dispatchable_demand, i), "q_demand_start"))

    # if bounded
    #     for (i, demand) in ref(wm, nw, :dispatchable_demand)
    #         flow_min, flow_max = 0.0, 0.0

    #         for name in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
    #             edges_fr = ref(wm, nw, Symbol(name * "_fr"), demand["node"])
    #             edges_to = ref(wm, nw, Symbol(name * "_to"), demand["node"])

    #             if length(edges_fr) > 0
    #                 flow_min -= sum(ref(wm, nw, Symbol(name), a)["flow_min"] for a in edges_fr)
    #                 flow_max -= sum(ref(wm, nw, Symbol(name), a)["flow_max"] for a in edges_fr)
    #             end

    #             if length(edges_to) > 0
    #                 flow_min -= sum(ref(wm, nw, Symbol(name), a)["flow_max"] for a in edges_to)
    #                 flow_max -= sum(ref(wm, nw, Symbol(name), a)["flow_min"] for a in edges_to)
    #             end
    #         end

    #         flow_min = max(flow_min, demand["flow_min"])
    #         flow_max = min(flow_max, demand["flow_max"])

    #         JuMP.set_lower_bound(q_demand[i], flow_min)
    #         JuMP.set_upper_bound(q_demand[i], flow_max)

    #         flow_mid = flow_min + 0.5 * (flow_max - flow_min)
    #         q_start = comp_start_value(demand, "q_demand_start", flow_mid)
    #         JuMP.set_start_value(q_demand[i], q_start)
    #     end
    # end

    report && sol_component_value(wm, nw, :demand, :q, ids(wm, nw, :dispatchable_demand), q_demand)
end


"""
    variable_tank_flow(
        wm::AbstractWaterModel;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true
    )

Creates bounded (by default) or unbounded volumetric flow rate variables
for tanks in the network at subnetwork (or time) index `nw`, i.e., `q_tank[i]`
for `i` in `tank`. Note that, unlike reservoirs, tanks can have inflow,
represented as a negative quantity.
"""
function variable_tank_flow(wm::AbstractWaterModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    q_tank = var(wm, nw)[:q_tank] = JuMP.@variable(wm.model,
        [i in ids(wm, nw, :tank)], base_name = "$(nw)_q_tank",
        start = comp_start_value(ref(wm, nw, :tank, i), "q_tank_start"))

    if bounded
        for (i, tank) in ref(wm, nw, :tank)
            flow_min, flow_max = 0.0, 0.0

            for name in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
                edges_fr = ref(wm, nw, Symbol(name * "_fr"), tank["node"])
                edges_to = ref(wm, nw, Symbol(name * "_to"), tank["node"])

                if length(edges_fr) > 0
                    flow_min += sum(ref(wm, nw, Symbol(name), a)["flow_min"] for a in edges_fr)
                    flow_max += sum(ref(wm, nw, Symbol(name), a)["flow_max"] for a in edges_fr)
                end

                if length(edges_to) > 0
                    flow_min -= sum(ref(wm, nw, Symbol(name), a)["flow_max"] for a in edges_to)
                    flow_max -= sum(ref(wm, nw, Symbol(name), a)["flow_min"] for a in edges_to)
                end
            end

            JuMP.set_lower_bound(q_tank[i], flow_min)
            JuMP.set_upper_bound(q_tank[i], flow_max)
            ref(wm, nw, :tank, i)["flow_min"] = flow_min
            ref(wm, nw, :tank, i)["flow_max"] = flow_max

            flow_mid = 0.5 * (flow_max + flow_min)
            q_start = comp_start_value(tank, "q_tank_start", flow_mid)
            JuMP.set_start_value(q_tank[i], q_start)
        end
    end

    report && sol_component_value(wm, nw, :tank, :q, ids(wm, nw, :tank), q_tank)
end


### Link variables. ###

"""
    variable_valve_indicator(
        wm::AbstractWaterModel;
        nw::Int=nw_id_default,
        relax::Bool=false,
        report::Bool=true
    )

Creates binary variables for valves in the network at subnetwork (or time)
index `nw`, i.e., `z_valve[a]` for `a` in `valve`. Here, one denotes that the
valve is open and zero denotes that the valve is closed.
"""
function variable_valve_indicator(wm::AbstractWaterModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if !relax
        z_valve = var(wm, nw)[:z_valve] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :valve)], base_name = "$(nw)_z_valve", binary = true,
            start = comp_start_value(ref(wm, nw, :valve, a), "z_valve_start", 1.0))
    else
        z_valve = var(wm, nw)[:z_valve] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :valve)], base_name = "z_valve[$(nw)]",
            lower_bound = 0.0, upper_bound = 1.0,
            start = comp_start_value(ref(wm, nw, :valve, a), "z_valve_start", 1.0))
    end

    for (a, valve) in ref(wm, nw, :valve)
        _fix_indicator_variable(z_valve[a], valve, "z")
    end

    report && sol_component_value(wm, nw, :valve, :status, ids(wm, nw, :valve), z_valve)
end


"""
    variable_regulator_indicator(
        wm::AbstractWaterModel;
        nw::Int=nw_id_default,
        relax::Bool=false,
        report::Bool=true
    )

Creates binary variables for regulators in the network at subnetwork (or time)
index `nw`, i.e., `z_regulator[a]` for `a` in `regulator`, where one denotes
that the pressure reducing regulator is currently active and zero otherwise.
"""
function variable_regulator_indicator(wm::AbstractWaterModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if !relax
        z_regulator = var(wm, nw)[:z_regulator] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :regulator)], base_name = "$(nw)_z_regulator", binary = true,
            start = comp_start_value(ref(wm, nw, :regulator, a), "z_regulator_start", 1.0))
    else
        z_regulator = var(wm, nw)[:z_regulator] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :regulator)], base_name = "$(nw)_z_regulator",
            lower_bound = 0.0, upper_bound = 1.0,
            start = comp_start_value(ref(wm, nw, :regulator, a), "z_regulator_start", 1.0))
    end

    for (a, regulator) in ref(wm, nw, :regulator)
        _fix_indicator_variable(z_regulator[a], regulator, "z")
    end

    report && sol_component_value(wm, nw, :regulator, :status, ids(wm, nw, :regulator), z_regulator)
end


"""
    variable_pump_indicator(
        wm::AbstractWaterModel;
        nw::Int=nw_id_default,
        relax::Bool=false,
        report::Bool=true
    )

Creates binary variables for pumps in the network at subnetwork (or time)
index `nw`, i.e., `z_pump[a]` for `a` in `pump`, where one denotes that the
pump is currently operating (i.e., on), and zero indicates that the pump is not
operating (i.e., off).
"""
function variable_pump_indicator(wm::AbstractWaterModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if !relax
        z_pump = var(wm, nw)[:z_pump] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pump)], base_name = "$(nw)_z_pump", binary = true,
            start = comp_start_value(ref(wm, nw, :pump, a), "z_pump_start", 1.0))
    else
        z_pump = var(wm, nw)[:z_pump] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pump)], base_name = "$(nw)_z_pump",
            lower_bound = 0.0, upper_bound = 1.0,
            start = comp_start_value(ref(wm, nw, :pump, a), "z_pump_start", 1.0))
    end

    for (a, pump) in ref(wm, nw, :pump)
        _fix_indicator_variable(z_pump[a], pump, "z")
    end

    report && sol_component_value(wm, nw, :pump, :status, ids(wm, nw, :pump), z_pump)
end


"""
    variable_pump_switch_on(
        wm::AbstractWaterModel;
        nw::Int=nw_id_default,
        relax::Bool=false,
        report::Bool=true
    )

Creates binary variables for pumps in the network at subnetwork (or time)
index `nw`, i.e., `z_switch_on_pump[a]` for `a` in `pump`, where one denotes
that the pump has been switched to the "on" status at the current subnetwork
(or time) index, and zero indicates that the pump has had no status change.
"""
function variable_pump_switch_on(wm::AbstractWaterModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if !relax
        z_switch_on_pump = var(wm, nw)[:z_switch_on_pump] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pump)], base_name = "$(nw)_z_switch_on_pump", binary = true,
            start = comp_start_value(ref(wm, nw, :pump, a), "z_switch_on_pump_start", 1.0))
    else
        z_switch_on_pump = var(wm, nw)[:z_switch_on_pump] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pump)], base_name = "$(nw)_z_switch_on_pump",
            lower_bound = 0.0, upper_bound = 1.0,
            start = comp_start_value(ref(wm, nw, :pump, a), "z_switch_on_pump_start", 1.0))
    end

    for (a, pump) in ref(wm, nw, :pump)
        _fix_indicator_variable(z_switch_on_pump[a], pump, "z_switch_on")
    end

    report && sol_component_value(wm, nw, :pump, :switch_on, ids(wm, nw, :pump), z_switch_on_pump)
end


"""
    variable_pump_switch_off(
        wm::AbstractWaterModel;
        nw::Int=nw_id_default,
        relax::Bool=false,
        report::Bool=true
    )

Creates binary variables for pumps in the network at subnetwork (or time)
index `nw`, i.e., `z_switch_off_pump[a]` for `a` in `pump`, where one denotes
that the pump has been switched to the "off" status at the current subnetwork
(or time) index, and zero indicates that the pump has had no status change.
"""
function variable_pump_switch_off(wm::AbstractWaterModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if !relax
        z_switch_off_pump = var(wm, nw)[:z_switch_off_pump] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pump)], base_name = "$(nw)_z_switch_off_pump", binary = true,
            start = comp_start_value(ref(wm, nw, :pump, a), "z_switch_off_pump_start", 1.0))
    else
        z_switch_off_pump = var(wm, nw)[:z_switch_off_pump] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pump)], base_name = "$(nw)_z_switch_off_pump",
            lower_bound = 0.0, upper_bound = 1.0,
            start = comp_start_value(ref(wm, nw, :pump, a), "z_switch_off_pump_start", 1.0))
    end

    for (a, pump) in ref(wm, nw, :pump)
        _fix_indicator_variable(z_switch_off_pump[a], pump, "z_switch_off")
    end

    report && sol_component_value(wm, nw, :pump, :switch_off, ids(wm, nw, :pump), z_switch_off_pump)
end


"""
    variable_des_pipe_indicator(
        wm::AbstractWaterModel;
        nw::Int=nw_id_default,
        relax::Bool=false,
        report::Bool=true
    )

Creates binary variables for design pipes in the network at subnetwork (or
time) index `nw`, i.e., `z_des_pipe[a]` for `a` in `des_pipe`, where one denotes
that the pipe has been selected within the design at the current subnetwork (or
time) index, and zero indicates that the design pipe has had not been selected.
"""
function variable_des_pipe_indicator(wm::AbstractWaterModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if !relax
        z_des_pipe = var(wm, nw)[:z_des_pipe] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :des_pipe)], base_name = "$(nw)_z_des_pipe", binary = true,
            start = comp_start_value(ref(wm, nw, :des_pipe, a), "z_des_pipe_start", 1.0))
    else
        z_des_pipe = var(wm, nw)[:z_des_pipe] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :des_pipe)], base_name = "$(nw)_z_des_pipe",
            lower_bound = 0.0, upper_bound = 1.0,
            start = comp_start_value(ref(wm, nw, :des_pipe, a), "z_des_pipe_start", 1.0))
    end

    for (a, des_pipe) in ref(wm, nw, :des_pipe)
        _fix_indicator_variable(z_des_pipe[a], des_pipe, "z")
    end

    report && sol_component_value(wm, nw, :des_pipe, :status, ids(wm, nw, :des_pipe), z_des_pipe)
end


function _fix_indicator_variable(v::JuMP.VariableRef, component::Dict{String, <:Any}, name::String)
    min_name, max_name = name * "_min", name * "_max"

    if haskey(component, min_name) && haskey(component, max_name)
        if component[min_name] == component[max_name]
            JuMP.fix(v, component[min_name])
        end
    end
end
