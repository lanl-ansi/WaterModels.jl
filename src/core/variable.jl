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
    h = var(wm, nw)[:h] = JuMP.@variable(wm.model,
        [i in ids(wm, nw, :node)], base_name="$(nw)_h",
        start=comp_start_value(ref(wm, nw, :node, i), "h_start"))
       
    if bounded
        h_lb, h_ub = calc_head_bounds(wm, nw)

        for (i, node) in ref(wm, nw, :node)
            JuMP.set_lower_bound(h[i], h_lb[i])
            JuMP.set_upper_bound(h[i], h_ub[i])
        end
    end

    report && sol_component_value(wm, nw, :node, :h, ids(wm, nw, :node), h)

    # Create expressions that calculate pressures.
    p = var(wm, nw)[:p] = JuMP.@expression(wm.model,
        [i in ids(wm, nw, :node)],
        var(wm, nw, :h, i) - ref(wm, nw, :node, i)["elevation"])

    report && sol_component_value(wm, nw, :node, :p, ids(wm, nw, :node), p)
end

"Creates head gain variables corresponding to all pumps in the network, i.e.,
`g[a]` for `a` in `pump`. These denote head gains between nodes i and j."
function variable_head_gain(wm::AbstractWaterModel; nw::Int=wm.cnw, report::Bool=true)
    g = var(wm, nw)[:g] = JuMP.@variable(wm.model, [a in ids(wm, nw, :pump)],
        base_name="$(nw)_g", lower_bound=0.0, # Pump flow is nonnegative.
        start=comp_start_value(ref(wm, nw, :pump, a), "g_start"))

    report && sol_component_value(wm, nw, :pump, :g, ids(wm, nw, :pump), g)
end

"Creates outgoing flow variables for all reservoirs in the network, i.e.,
`qr[i]` for `i` in `reservoir`. Note that these variables are always
nonnegative, as there is never incoming flow to a reservoir."
function variable_reservoir(wm::AbstractWaterModel; nw::Int=wm.cnw, report::Bool=true)
    qr = var(wm, nw)[:qr] = JuMP.@variable(wm.model,
        [i in ids(wm, nw, :reservoir)], lower_bound=0.0, base_name="$(nw)_qr",
        start=comp_start_value(ref(wm, nw, :reservoir, i), "qr_start"))

    report && sol_component_value(wm, nw, :reservoir, :qr, ids(wm, nw, :reservoir), qr)
end

"Creates outgoing flow variables for all tanks in the network, i.e., `qt[i]`
for `i` in `tank`. Note that, unlike reservoirs, tanks can have inflow."
function variable_tank(wm::AbstractWaterModel; nw::Int=wm.cnw, report::Bool=true)
    qt = var(wm, nw)[:qt] = JuMP.@variable(wm.model,
        [i in ids(wm, nw, :tank)], base_name="$(nw)_qt",
        start=comp_start_value(ref(wm, nw, :tank, i), "qt_start"))

    report && sol_component_value(wm, nw, :tank, :qt, ids(wm, nw, :tank), qt)

    V = var(wm, nw)[:V] = JuMP.@expression(wm.model,
        [i in ids(wm, nw, :tank)],
        (var(wm, nw, :h, ref(wm, nw, :tank, i)["node"]) -
         ref(wm, nw, :node, ref(wm, nw, :tank, i)["node"])["elevation"])
         * 0.25 * pi * ref(wm, nw, :tank, i)["diameter"]^2)

    report && sol_component_value(wm, nw, :tank, :V, ids(wm, nw, :tank), V)
end


### Link variables. ###
"Creates binary variables for all check valves in the network, i.e.,
`z_check_valve[a]` for `a` in `check_valve`, where one denotes that the check
valve is open and zero denotes that the check valve is closed."
function variable_check_valve_indicator(wm::AbstractWaterModel; nw::Int=wm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_check_valve = var(wm, nw)[:z_check_valve] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :check_valve)], base_name = "$(nw)_z_check_valve",
            binary = true,
            start = comp_start_value(ref(wm, nw, :check_valve, a), "z_check_valve_start"))
    else
        z_check_valve = var(wm, nw)[:z_check_valve] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :check_valve)], base_name = "z_check_valve[$(nw)]",
            lower_bound = 0.0, upper_bound = 1.0,
            start = comp_start_value(ref(wm, nw, :check_valve, a), "z_check_valve_start"))
    end

    report && _IM.sol_component_value(wm, nw, :check_valve, :status, ids(wm, nw, :check_valve), z_check_valve)
end

"Creates binary variables for all shutoff valves in the network, i.e.,
`z_shutoff_valve[a]` for `a` in `sv`, where one denotes that the shutoff valve
is open and zero denotes that the shutoff valve is closed."
function variable_shutoff_valve_indicator(wm::AbstractWaterModel; nw::Int=wm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_shutoff_valve = var(wm, nw)[:z_shutoff_valve] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :shutoff_valve)], base_name = "$(nw)_z_shutoff_valve",
            binary = true,
            start = comp_start_value(ref(wm, nw, :shutoff_valve, a), "z_shutoff_valve_start"))
    else
        z_shutoff_valve = var(wm, nw)[:z_shutoff_valve] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :shutoff_valve)], base_name = "$(nw)_z_shutoff_valve",
            lower_bound = 0.0, upper_bound = 1.0,
            start = comp_start_value(ref(wm, nw, :shutoff_valve, a), "z_shutoff_valve_start"))
    end

    report && sol_component_value(wm, nw, :shutoff_valve, :status, ids(wm, nw, :shutoff_valve), z_shutoff_valve)
end

"Creates binary variables for all PRVs in the network, i.e.,
`z_pressure_reducing_valve[a]` for `a` in `pressure_reducing_valve`, where one
denotes that the pressure reducing is currently open and zero otherwise."
function variable_pressure_reducing_valve_indicator(wm::AbstractWaterModel; nw::Int=wm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_pressure_reducing_valve = var(wm, nw)[:z_pressure_reducing_valve] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pressure_reducing_valve)], base_name =
            "$(nw)_z_pressure_reducing_valve", binary = true, start =
            comp_start_value(ref(wm, nw, :pressure_reducing_valve, a),
                             "z_pressure_reducing_valve_start"))
    else
        z_pressure_reducing_valve = var(wm, nw)[:z_pressure_reducing_valve] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pressure_reducing_valve)], base_name = "$(nw)_z_pressure_reducing_valve",
            lower_bound = 0.0, upper_bound = 1.0,
            start = comp_start_value(ref(wm, nw, :pressure_reducing_valve, a), "z_pressure_reducing_valve_start"))
    end

    report && sol_component_value(wm, nw, :pressure_reducing_valve, :z_pressure_reducing_valve,
        ids(wm, nw, :pressure_reducing_valve), z_pressure_reducing_valve)
end

"Creates binary variables for all pumps in the network, i.e., `z_pump[a]` for
`a` in `pump`, where one denotes that the pump is currently operating (i.e.,
on), and zero indicates that the pump is not operating (i.e., off)."
function variable_pump_indicator(wm::AbstractWaterModel; nw::Int=wm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_pump = var(wm, nw)[:z_pump] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pump)], base_name = "$(nw)_z_pump",
            binary = true,
            start = comp_start_value(ref(wm, nw, :pump, a), "z_pump_start"))
    else
        z_pump = var(wm, nw)[:z_pump] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pump)], base_name = "$(nw)_z_pump",
            lower_bound = 0.0, upper_bound = 1.0,
            start = comp_start_value(ref(wm, nw, :pump, a), "z_pump_start"))
    end

    report && sol_component_value(wm, nw, :pump, :status, ids(wm, nw, :pump), z_pump)
end

function variable_pressure_reducing_valve_operation(wm::AbstractWaterModel; nw::Int=wm.cnw, report::Bool=true)
    variable_pressure_reducing_valve_indicator(wm, nw=nw, report=report)
end

function variable_pump_operation(wm::AbstractWaterModel; nw::Int=wm.cnw, report::Bool=true)
end

"Creates binary variables for all network design or design resistances in
the network, i.e., `x_res[a]` for `a` in `pipe`, for `r` in `resistance[a]`,
where one denotes that the given resistance is active in the design."
function variable_resistance(wm::AbstractWaterModel; nw::Int=wm.cnw, report::Bool=true)
    x_res = var(wm, nw)[:x_res] = Dict{Int, Array{JuMP.VariableRef}}()

    for a in ids(wm, nw, :pipe_des)
        n_r = length(ref(wm, nw, :resistance, a)) # Number of resistances.
        var(wm, nw, :x_res)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            binary=true, base_name="$(nw)_x_res[$(a)]",
            start=comp_start_value(ref(wm, nw, :pipe_des, a), "x_res_start", r))
    end

    report && sol_component_value(wm, nw, :pipe_des, :x_res, ids(wm, nw, :pipe_des), x_res)
end
