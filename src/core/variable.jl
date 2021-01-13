########################################################################
# This file defines commonly-used variables for water systems models.
########################################################################


"Sets the start value for a given variable."
function comp_start_value(comp::Dict{String,<:Any}, key::String, default::Float64=0.0)
    return get(comp, key, default)
end


"Sets the start value for a given variable."
function comp_start_value(comp::Dict{String,<:Any}, key_1::String, key_2::Int64, default=0.0)
    return key_1 in keys(comp) ? get(get(comp, key_1, default), key_2, default) : default
end


"Given a variable that is indexed by component IDs, builds the standard solution structure."
function sol_component_value(wm::AbstractWaterModel, n::Int, comp_name::Symbol, field_name::Symbol, comp_ids, variables)
    for i in comp_ids
        @assert !haskey(sol(wm, n, comp_name, i), field_name)
        sol(wm, n, comp_name, i)[field_name] = variables[i]
    end
end


### Variables related to nodal components. ###
"Creates bounded (by default) or unbounded total hydraulic head (or head)
variables for all nodes in the network, i.e., `h[i]` for `i` in `node`."
function variable_head(wm::AbstractWaterModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    # Initialize variables for total hydraulic head.
    h = var(wm, nw)[:h] = JuMP.@variable(wm.model,
        [i in ids(wm, nw, :node)], base_name = "$(nw)_h",
        start = comp_start_value(ref(wm, nw, :node, i), "h_start"))
       
    if bounded
        for (i, node) in ref(wm, nw, :node)
            # Set the lower and upper bounds for each head.
            JuMP.set_lower_bound(h[i], node["head_min"])
            JuMP.set_upper_bound(h[i], node["head_max"])

            # Set start value for the head variable with possibly better data.
            h_mid = node["head_min"] + 0.5 * (node["head_max"] - node["head_min"])
            h_start = comp_start_value(node, "h_start", h_mid)
            JuMP.set_start_value(h[i], h_start)
        end
    end

    # Initialize an entry to the solution component dictionary for heads.
    report && sol_component_value(wm, nw, :node, :h, ids(wm, nw, :node), h)

    # Create expressions that calculate pressures from head.
    p = var(wm, nw)[:p] = JuMP.@expression(wm.model,
        [i in ids(wm, nw, :node)], var(wm, nw, :h, i) - ref(wm, nw, :node, i)["elevation"])

    # Initialize an entry to the solution component dictionary for pressures.
    report && sol_component_value(wm, nw, :node, :p, ids(wm, nw, :node), p)

    # Create an expression that maps total hydraulic head to volume for tanks.
    V = var(wm, nw)[:V] = JuMP.@expression(wm.model, [i in ids(wm, nw, :tank)],
        (var(wm, nw, :h, ref(wm, nw, :tank, i)["node"]) -
         ref(wm, nw, :node, ref(wm, nw, :tank, i)["node"])["elevation"])
         * 0.25 * pi * ref(wm, nw, :tank, i)["diameter"]^2)

    # Initialize an entry to the solution component dictionary for volumes.
    report && sol_component_value(wm, nw, :tank, :V, ids(wm, nw, :tank), V)
end


function variable_pump_head_gain(wm::AbstractWaterModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
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


function variable_pump_power(wm::AbstractWaterModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    # Initialize variables for the power utilization of a pump.
    Ps = var(wm, nw)[:Ps_pump] = JuMP.@variable(wm.model, [a in ids(wm, nw, :pump)],
        base_name="$(nw)_Ps_pump", lower_bound=0.0, # Pump power is nonnegative.
        start=comp_start_value(ref(wm, nw, :pump, a), "Ps_pump_start", 0.0))

    if bounded
        for (a, pump) in ref(wm, nw, :pump)
            # Get the nodes that are connected by the pump.
            node_fr = ref(wm, nw, :node, pump["node_fr"])
            node_to = ref(wm, nw, :node, pump["node_to"])
    
            # Set the upper bound for the power variable.
            P_max = _calc_pump_power_max(pump, node_fr, node_to)
            JuMP.set_upper_bound(Ps[a], P_max / (_DENSITY * _GRAVITY))

            # Set the start value for the scaled power variable.
            Ps_mid = 0.5 * P_max / (_DENSITY * _GRAVITY)
            Ps_start = comp_start_value(pump, "Ps_pump_start", Ps_mid)
            JuMP.set_start_value(Ps[a], Ps_start)
        end
    end

    # Initialize an entry to the solution component dictionary for powers.
    report && sol_component_value(wm, nw, :pump, :Ps, ids(wm, nw, :pump), Ps)

    # Create expressions to compute the unscaled power values.
    P = var(wm, nw)[:P_pump] = JuMP.@expression(
        wm.model, [a in ids(wm, nw, :pump)], Ps[a] * _DENSITY * _GRAVITY)

    # Initialize an entry to the solution component dictionary for powers.
    report && sol_component_value(wm, nw, :pump, :P, ids(wm, nw, :pump), P)
end


"Instantiates outgoing flow variables for all reservoirs in the network, i.e.,
`q_reservoir[i]` for `i` in `reservoir`. Note that these variables are always nonnegative,
since for each reservoir, there will never be incoming flow."
function variable_reservoir_flow(wm::AbstractWaterModel; nw::Int=wm.cnw, report::Bool=true)
    q_reservoir = var(wm, nw)[:q_reservoir] = JuMP.@variable(wm.model,
        [i in ids(wm, nw, :reservoir)], lower_bound=0.0, base_name="$(nw)_q_reservoir",
        start=comp_start_value(ref(wm, nw, :reservoir, i), "q_reservoir_start"))

    report && sol_component_value(wm, nw, :reservoir, :q, ids(wm, nw, :reservoir), q_reservoir)
end


"Instantiates demand flow variables for all dispatchable demands in the network, i.e.,
`demand[i]` for `i` in `dispatchable_demand`."
function variable_demand_flow(wm::AbstractWaterModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    q_demand = var(wm, nw)[:q_demand] = JuMP.@variable(wm.model,
        [i in ids(wm, nw, :dispatchable_demand)], base_name="$(nw)_q_demand",
        start=comp_start_value(ref(wm, nw, :dispatchable_demand, i), "q_demand_start"))

    if bounded
        for (i, demand) in ref(wm, nw, :dispatchable_demand)
            JuMP.set_lower_bound(q_demand[i], demand["flow_min"])
            JuMP.set_upper_bound(q_demand[i], demand["flow_max"])

            # Set potentially better start value based on demand bounds.
            q_mid = demand["flow_min"] + 0.5 * (demand["flow_max"] - demand["flow_min"])
            q_start = comp_start_value(demand, "q_demand_start", q_mid)
            JuMP.set_start_value(q_demand[i], q_start)
        end
    end

    report && sol_component_value(wm, nw, :demand, :q, ids(wm, nw, :dispatchable_demand), q_demand)
end


"Creates outgoing flow variables for all tanks in the network, i.e., `q_tank[i]`
for `i` in `tank`. Note that, unlike reservoirs, tanks can have inflow."
function variable_tank_flow(wm::AbstractWaterModel; nw::Int=wm.cnw, report::Bool=true)
    q_tank = var(wm, nw)[:q_tank] = JuMP.@variable(wm.model,
        [i in ids(wm, nw, :tank)], base_name="$(nw)_q_tank",
        start=comp_start_value(ref(wm, nw, :tank, i), "q_tank_start"))

    report && sol_component_value(wm, nw, :tank, :q, ids(wm, nw, :tank), q_tank)
end


### Link variables. ###
"Creates binary variables for valves in the network, i.e., `z_valve[a]` for `a` in `valve`,
where one denotes that the valve is open and zero denotes that the valve is closed."
function variable_valve_indicator(wm::AbstractWaterModel; nw::Int=wm.cnw, relax::Bool=false, report::Bool=true)
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

    report && _IM.sol_component_value(wm, nw, :valve, :status, ids(wm, nw, :valve), z_valve)
end


"Creates binary variables for all regulators in the network, i.e., `z_regulator[a]` for `a` in
`regulator`, where one denotes that the pressure reducing is currently open and zero otherwise."
function variable_regulator_indicator(wm::AbstractWaterModel; nw::Int=wm.cnw, relax::Bool=false, report::Bool=true)
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


"Creates binary variables for all pumps in the network, i.e., `z_pump[a]` for
`a` in `pump`, where one denotes that the pump is currently operating (i.e.,
on), and zero indicates that the pump is not operating (i.e., off)."
function variable_pump_indicator(wm::AbstractWaterModel; nw::Int=wm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_pump = var(wm, nw)[:z_pump] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pump)], base_name = "$(nw)_z_pump",
            binary = true,
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


"Creates binary variables for all design pipes in the network, i.e.,
`z_des_pipe[a]` for `a` in `des_pipe`, where one denotes that the pipe is
selected within the design, and zero denotes that the pipe is not selected."
function variable_des_pipe_indicator(wm::AbstractWaterModel; nw::Int=wm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_des_pipe = var(wm, nw)[:z_des_pipe] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :des_pipe)], base_name = "$(nw)_z_des_pipe",
            binary = true,
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
