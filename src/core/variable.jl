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

"Given a constant value, builds the standard componentwise solution structure."
function sol_component_fixed(wm::AbstractWaterModel, n::Int, comp_name::Symbol, field_name::Symbol, comp_ids, constant)
    for i in comp_ids
        @assert !haskey(sol(wm, n, comp_name, i), field_name)
        sol(wm, n, comp_name, i)[field_name] = constant
    end
end

"Given a variable that is indexed by component IDs, builds the standard solution structure."
function sol_component_value(wm::AbstractWaterModel, n::Int, comp_name::Symbol, field_name::Symbol, comp_ids, variables)
    for i in comp_ids
        @assert !haskey(sol(wm, n, comp_name, i), field_name)
        sol(wm, n, comp_name, i)[field_name] = variables[i]
    end
end

"Maps asymmetric edge variables into components."
function sol_component_value_edge(wm::AbstractWaterModel, n::Int, comp_name::Symbol, field_name_fr::Symbol, field_name_to::Symbol, comp_ids_fr, comp_ids_to, variables)
    for (l, i, j) in comp_ids_fr
        @assert !haskey(sol(wm, n, comp_name, l), field_name_fr)
        sol(wm, n, comp_name, l)[field_name_fr] = variables[(l, i, j)]
    end

    for (l, i, j) in comp_ids_to
        @assert !haskey(sol(wm, n, comp_name, l), field_name_to)
        sol(wm, n, comp_name, l)[field_name_to] = variables[(l, i, j)]
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
    qt = var(wm, nw)[:qt] = JuMP.@variable(wm.model, [i in ids(wm, nw, :tank)],
        base_name="$(nw)_qt",
        start=comp_start_value(ref(wm, nw, :tank, i), "qt_start"))

    report && sol_component_value(wm, nw, :tank, :qt, ids(wm, nw, :tank), qt)
end

"Creates bounded (by default) or unbounded (but still nonnegative) volume
variables for all tanks in the network, i.e., `V[i]` for `i` in `tank`."
function variable_volume(wm::AbstractWaterModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    V = var(wm, nw)[:V] = JuMP.@variable(wm.model,
        [i in ids(wm, nw, :tank)], base_name="$(nw)_V", lower_bound=0.0,
        start=comp_start_value(ref(wm, nw, :tank, i), "V_start"))
       
    if bounded
        V_lb, V_ub = calc_tank_volume_bounds(wm, nw)

        for (i, tank) in ref(wm, nw, :tank)
            JuMP.set_lower_bound(V[i], V_lb[i])
            JuMP.set_upper_bound(V[i], V_ub[i])
        end
    end

    report && sol_component_value(wm, nw, :tank, :V, ids(wm, nw, :tank), V)
end

### Link variables. ###
"Creates binary variables for all check valves in the network, i.e., `x_cv[a]`
for `a` in `check_valve`, where one denotes that the check valve is open."
function variable_check_valve(wm::AbstractWaterModel; nw::Int=wm.cnw, report::Bool=true)
    x_cv = var(wm, nw)[:x_cv] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :check_valve)], base_name="x_cv[$(nw)]", binary=true,
        start=comp_start_value(ref(wm, nw, :check_valve, a), "x_cv_start"))

    report && sol_component_value(wm, nw, :check_valve, :x_cv,
        ids(wm, nw, :check_valve), x_cv)
end

"Creates binary variables for all pumps in the network, i.e., `x_pump[a]`
for `a` in `pump`, where one denotes that the pump is currently on."
function variable_fixed_speed_pump_operation(wm::AbstractWaterModel; nw::Int=wm.cnw, report::Bool=true)
    x_pump = var(wm, nw)[:x_pump] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pump)], base_name="x_pump[$(nw)]", binary=true,
        start=comp_start_value(ref(wm, nw, :pump, a), "x_pump_start"))

    report && sol_component_value(wm, nw, :pump, :x_pump, ids(wm, nw, :pump), x_pump)
end

"Creates binary variables for all pumps in the network, i.e., `x_pump[a]`
for `a` in `pump`, where one denotes that the pump is currently on."
function variable_fixed_speed_pump_threshold(wm::AbstractWaterModel; nw::Int=wm.cnw)
    var(wm, nw)[:x_thrs_gt] = JuMP.@variable(wm.model, [a in ids(wm, nw, :pump)],
        base_name="x_thrs_gt[$(nw)]", binary=true,
        start=comp_start_value(ref(wm, nw, :pump, a), "x_thrs_gt_start"))
    var(wm, nw)[:x_thrs_lt] = JuMP.@variable(wm.model, [a in ids(wm, nw, :pump)],
        base_name="x_thrs_lt[$(nw)]", binary=true,
        start=comp_start_value(ref(wm, nw, :pump, a), "x_thrs_lt_start"))
    var(wm, nw)[:x_thrs_bt] = JuMP.@variable(wm.model, [a in ids(wm, nw, :pump)],
        base_name="x_thrs_bt[$(nw)]", binary=true,
        start=comp_start_value(ref(wm, nw, :pump, a), "x_thrs_bt_start"))
end

function variable_pump_operation(wm::AbstractWaterModel; nw::Int=wm.cnw)
    variable_fixed_speed_pump_operation(wm, nw=nw)
end

function variable_pump_control(wm::AbstractWaterModel; nw::Int=wm.cnw)
    variable_fixed_speed_pump_threshold(wm, nw=nw)
end

"Creates binary variables for all network expansion or design resistances in
the network, i.e., `x_res[a]` for `a` in `pipe`, for `r` in `resistance[a]`,
where one denotes that the given resistance is active in the design."
function variable_resistance(wm::AbstractWaterModel; nw::Int=wm.cnw, report::Bool=true)
    x_res = var(wm, nw)[:x_res] = Dict{Int, Array{JuMP.VariableRef}}()

    for a in ids(wm, nw, :link_ne)
        n_r = length(ref(wm, nw, :resistance, a)) # Number of resistances.
        var(wm, nw, :x_res)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            binary=true, base_name="$(nw)_x_res[$(a)]",
            start=comp_start_value(ref(wm, nw, :link_ne, a), "x_res_start", r))
    end

    report && sol_component_value(wm, nw, :pipe_ne, :x_res, ids(wm, nw, :pipe_ne), x_res)
end
