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
        [i in ids(wm, nw, :node)], base_name="$(nw)_h",
        start=comp_start_value(ref(wm, nw, :node, i), "h_start"))
       
    if bounded
        # Compute the lower and upper head bound dictionaries.
        h_lb, h_ub = calc_head_bounds(wm, nw)

        for (i, node) in ref(wm, nw, :node)
            # Set the lower and upper bounds for each head.
            JuMP.set_lower_bound(h[i], h_lb[i])
            JuMP.set_upper_bound(h[i], h_ub[i])
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

    # Initialize an entry to the solution component dictionary for head gains.
    report && sol_component_value(wm, nw, :pump, :g, ids(wm, nw, :pump), g)
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
        q_demand_lb, q_demand_ub = calc_demand_bounds(wm, nw)

        for (i, demand) in ref(wm, nw, :dispatchable_demand)
            JuMP.set_lower_bound(q_demand[i], q_demand_lb[i])
            JuMP.set_upper_bound(q_demand[i], q_demand_ub[i])
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
            start = comp_start_value(ref(wm, nw, :valve, a), "z_valve_start"))
    else
        z_valve = var(wm, nw)[:z_valve] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :valve)], base_name = "z_valve[$(nw)]",
            lower_bound = 0.0, upper_bound = 1.0,
            start = comp_start_value(ref(wm, nw, :valve, a), "z_valve_start"))
    end

    report && _IM.sol_component_value(wm, nw, :valve, :status, ids(wm, nw, :valve), z_valve)
end


"Creates binary variables for all regulators in the network, i.e., `z_regulator[a]` for `a` in
`regulator`, where one denotes that the pressure reducing is currently open and zero otherwise."
function variable_regulator_indicator(wm::AbstractWaterModel; nw::Int=wm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_regulator = var(wm, nw)[:z_regulator] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :regulator)], base_name = "$(nw)_z_regulator", binary = true,
            start = comp_start_value(ref(wm, nw, :regulator, a), "z_regulator_start"))
    else
        z_regulator = var(wm, nw)[:z_regulator] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :regulator)], base_name = "$(nw)_z_regulator",
            lower_bound = 0.0, upper_bound = 1.0,
            start = comp_start_value(ref(wm, nw, :regulator, a), "z_regulator_start"))
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
            binary = true, start = comp_start_value(ref(wm, nw, :pump, a), "z_pump_start"))
    else
        z_pump = var(wm, nw)[:z_pump] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pump)], base_name = "$(nw)_z_pump",
            lower_bound = 0.0, upper_bound = 1.0,
            start = comp_start_value(ref(wm, nw, :pump, a), "z_pump_start"))
    end

    report && sol_component_value(wm, nw, :pump, :status, ids(wm, nw, :pump), z_pump)
end


"Creates binary variables for all network design or design resistances in
the network, i.e., `x_res[a]` for `a` in `pipe`, for `r` in `resistance[a]`,
where one denotes that the given resistance is active in the design."
function variable_resistance(wm::AbstractWaterModel; nw::Int=wm.cnw, report::Bool=true)
    x_res = var(wm, nw)[:x_res] = Dict{Int, Array{JuMP.VariableRef}}()

    for a in ids(wm, nw, :des_pipe)
        n_r = length(ref(wm, nw, :resistance, a)) # Number of resistances.
        var(wm, nw, :x_res)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            binary=true, base_name="$(nw)_x_res[$(a)]",
            start=comp_start_value(ref(wm, nw, :des_pipe, a), "x_res_start", r))
    end

    report && sol_component_value(wm, nw, :des_pipe, :x_res, ids(wm, nw, :des_pipe), x_res)
end
