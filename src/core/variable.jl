########################################################################
# This file defines commonly-used variables for water systems models.
########################################################################

function get_start(set, item_key, value_key, default)
    return get(get(set, item_key, Dict()), value_key, default)
end

function get_start(set, item_key, first_value_key, second_value_key, default)
    return get(get(get(set, item_key, Dict()), second_value_key, Dict()), first_value_key, default)
end

### Variables related to nodal components. ###
"Creates bounded (by default) or unbounded total hydraulic head variables for
all nodes in the network, i.e., `h[i]` for `i` in `node`."
function variable_head(wm::AbstractWaterModel; nw::Int=wm.cnw, bounded::Bool=true)
    bounded ? variable_head_bounded(wm, nw=nw) :
        variable_head_unbounded(wm, nw=nw)
end

"Creates bounded total hydraulic head variables for all nodes in the network,
i.e., `h[i]` for `i` in `node`."
function variable_head_bounded(wm::AbstractWaterModel; nw::Int=wm.cnw)
    h_lb, h_ub = calc_head_bounds(wm, nw)
    var(wm, nw)[:h] = JuMP.@variable(wm.model, [i in ids(wm, nw, :node)],
        base_name="h[$(nw)]", lower_bound=h_lb[i], upper_bound=h_ub[i],
        start=get_start(ref(wm, nw, :node), i, "h_start", h_ub[i]))
end

"Creates head gain variables corresponding to all pumps in the network, i.e.,
`g[a]` for `a` in `pump`. These denote head gains between nodes i and j."
function variable_head_gain(wm::AbstractWaterModel; nw::Int=wm.cnw)
    var(wm, nw)[:g] = JuMP.@variable(wm.model, [a in ids(wm, nw, :pump)],
        base_name="g[$(nw)]", lower_bound=0.0,
        start=get_start(ref(wm, nw, :pump), a, "g", 0.0))
end

"Creates unbounded total hydraulic head variables for all nodes in the network,
i.e., `h[i]` for `i` in `node`. Note that head variables corresponding to
reservoir and tank component types are also unbounded, here."
function variable_head_unbounded(wm::AbstractWaterModel; nw::Int=wm.cnw)
    var(wm, nw)[:h] = JuMP.@variable(wm.model, [i in ids(wm, nw, :node)],
        base_name="h[$(nw)]",
        start=get_start(ref(wm, nw, :node), i, "h_start", 0.0))
end

"Creates variables denoting absolute constraint violations with respect to head
bounds for all nodes in the network, i.e., `h_viol[i]` for `i` in `node`."
function variable_head_violation(wm::AbstractWaterModel; nw::Int=wm.cnw)
    var(wm, nw)[:h_viol] = JuMP.@variable(wm.model, [i in ids(wm, nw, :node)],
        base_name="h_viol[$(nw)]", lower_bound=0.0,
        start=get_start(ref(wm, nw, :node), i, "h_viol_start", 0.0))
end

"Creates outgoing flow variables for all reservoirs in the network, i.e.,
`qr[i]` for `i` in `reservoir`. Note that these variables are always
nonnegative, as there is never incoming flow to a reservoir."
function variable_reservoir(wm::AbstractWaterModel; nw::Int=wm.cnw)
    var(wm, nw)[:qr] = JuMP.@variable(wm.model, [i in ids(wm, nw, :reservoir)],
        base_name="qr[$(nw)]", lower_bound=0.0,
        start=get_start(ref(wm, nw, :reservoir), i, "qr_start", 0.0))
end

"Creates outgoing flow variables for all tanks in the network, i.e., `qt[i]`
for `i` in `tank`. Note that, unlike reservoirs, tanks can have inflow."
function variable_tank(wm::AbstractWaterModel; nw::Int=wm.cnw)
    var(wm, nw)[:qt] = JuMP.@variable(wm.model, [i in ids(wm, nw, :tank)],
        base_name="qt[$(nw)]",
        start=get_start(ref(wm, nw, :node), i, "qt_start", 0.0))
end

"Creates bounded (by default) or unbounded (but still nonnegative) volume
variables for all tanks in the network, i.e., `V[i]` for `i` in `tank`."
function variable_volume(wm::AbstractWaterModel; nw::Int=wm.cnw, bounded::Bool=true)
    bounded ? variable_volume_bounded(wm, nw=nw) :
        variable_volume_unbounded(wm, nw=nw)
end

"Creates bounded volume variables for all tanks in the network, i.e., `V[i]`
for `i` in `tank`. Bounds are related to minimum and maximum tank levels."
function variable_volume_bounded(wm::AbstractWaterModel; nw::Int=wm.cnw)
    V_lb, V_ub = calc_tank_volume_bounds(wm, nw)
    var(wm, nw)[:V] = JuMP.@variable(wm.model, [i in ids(wm, nw, :tank)],
        base_name="V[$(nw)]", lower_bound=V_lb[i], upper_bound=V_ub[i],
        start=get_start(ref(wm, nw, :tank), i, "V_start", V_ub[i]))
end

"Creates unbounded volume variables for all tanks in the network, i.e., `V[i]`
for `i` in `tank`. Note that since volume is inherently a nonnegative quantity,
however, these volume variables still contain lower bounds of zero."
function variable_volume_unbounded(wm::AbstractWaterModel; nw::Int=wm.cnw)
    var(wm, nw)[:V] = JuMP.@variable(wm.model, [i in ids(wm, nw, :tank)],
        base_name="V[$(nw)]", lower_bound=0.0,
        start=get_start(ref(wm, nw, :tank), i, "V_start", 0.0))
end

### Link variables. ###
function variable_check_valve(wm::AbstractWaterModel; nw::Int=wm.cnw)
    # Create variables for the statuses of check valves (one is open).
    var(wm, nw)[:x_cv] = JuMP.@variable(wm.model,
        [i in ids(wm, nw, :check_valve)], base_name="x_cv[$(nw)]", binary=true,
        start=get_start(ref(wm, nw, :link), i, "x_cv_start", 0))
end

function variable_flow_violation(wm::AbstractWaterModel; nw::Int=wm.cnw)
    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    var(wm, nw)[:Delta_e] = JuMP.@variable(wm.model, [a in ids(wm, nw, :link)],
        base_name="Delta_e[$(nw)]", lower_bound=0.0,
        start=get_start(ref(wm, nw, :link), a, "Delta_e_start", 0.0))
end

function variable_fixed_speed_pump_operation(wm::AbstractWaterModel; nw::Int=wm.cnw)
    var(wm, nw)[:x_pump] = JuMP.@variable(wm.model, [a in ids(wm, nw, :pump)],
        base_name="x_pump[$(nw)]", binary=true,
        start=get_start(ref(wm, nw, :pump), a, "x_pump_start", 1.0))
end

function variable_fixed_speed_pump_threshold(wm::AbstractWaterModel; nw::Int=wm.cnw)
    var(wm, nw)[:x_thrs_gt] = JuMP.@variable(wm.model, [a in ids(wm, nw, :pump)],
        base_name="x_thrs_gt[$(nw)]", binary=true,
        start=get_start(ref(wm, nw, :pump), a, "x_thrs_gt_start", 0.0))
    var(wm, nw)[:x_thrs_lt] = JuMP.@variable(wm.model, [a in ids(wm, nw, :pump)],
        base_name="x_thrs_lt[$(nw)]", binary=true,
        start=get_start(ref(wm, nw, :pump), a, "x_thrs_lt_start", 0.0))
    var(wm, nw)[:x_thrs_bt] = JuMP.@variable(wm.model, [a in ids(wm, nw, :pump)],
        base_name="x_thrs_bt[$(nw)]", binary=true,
        start=get_start(ref(wm, nw, :pump), a, "x_thrs_bt_start", 1.0))
end

function variable_pump_operation(wm::AbstractWaterModel; nw::Int=wm.cnw)
    variable_fixed_speed_pump_operation(wm, nw=nw)
end

function variable_pump_control(wm::AbstractWaterModel; nw::Int=wm.cnw)
    variable_fixed_speed_pump_threshold(wm, nw=nw)
end

function variable_resistance(wm::AbstractWaterModel; nw::Int=wm.cnw)
    link_ids = ids(wm, nw, :link_ne)
    var(wm, nw)[:x_res] = Dict{Int, Array{JuMP.VariableRef}}(a => [] for a in link_ids)

    for a in ids(wm, nw, :link_ne)
        n_r = length(ref(wm, nw, :resistance, a)) # Number of resistances.
        var(wm, nw, :x_res)[a] = JuMP.@variable(wm.model, [r in 1:n_r],
            binary=true, base_name="x_res[$(nw)][$(a)]",
            start=get_start(ref(wm, nw, :link_ne), a, r, "x_res_start", 0.0))
    end
end
