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
function variable_head(wm::AbstractWaterModel, n::Int=wm.cnw; bounded::Bool=true)
    bounded ? variable_head_bounded(wm, n) : variable_head_unbounded(wm, n)
end

"Creates bounded total hydraulic head variables for all nodes in the network,
i.e., `h[i]` for `i` in `node`."
function variable_head_bounded(wm::AbstractWaterModel, n::Int=wm.cnw)
    h_lb, h_ub = calc_head_bounds(wm, n)
    var(wm, n)[:h] = JuMP.@variable(wm.model, [i in ids(wm, n, :node)],
        base_name="h[$(n)]", lower_bound=h_lb[i], upper_bound=h_ub[i],
        start=get_start(ref(wm, n, :node), i, "h_start", h_ub[i]))
end

"Creates head gain variables corresponding to all pumps in the network, i.e.,
`g[a]` for `a` in `pump`. These denote head gains between nodes i and j."
function variable_head_gain(wm::AbstractWaterModel, n::Int=wm.cnw)
    var(wm, n)[:g] = JuMP.@variable(wm.model, [a in ids(wm, n, :pump)],
        base_name="g[$(n)]", lower_bound=0.0,
        start=get_start(ref(wm, n, :pump), a, "g", 0.0))
end

"Creates unbounded total hydraulic head variables for all nodes in the network,
i.e., `h[i]` for `i` in `node`. Note that head variables corresponding to
reservoir and tank component types are also unbounded, here."
function variable_head_unbounded(wm::AbstractWaterModel, n::Int=wm.cnw)
    var(wm, n)[:h] = JuMP.@variable(wm.model, [i in ids(wm, n, :node)],
        base_name="h[$(n)]",
        start=get_start(ref(wm, n, :node), i, "h_start", 0.0))
end

"Creates variables denoting absolute constraint violations with respect to head
bounds for all nodes in the network, i.e., `h_viol[i]` for `i` in `node`."
function variable_head_violation(wm::AbstractWaterModel, n::Int=wm.cnw)
    var(wm, n)[:h_viol] = JuMP.@variable(wm.model, [i in ids(wm, n, :node)],
        base_name="h_viol[$(n)]", lower_bound=0.0,
        start=get_start(ref(wm, n, :node), i, "h_viol_start", 0.0))
end

"Creates outgoing flow variables for all reservoirs in the network, i.e.,
`qr[i]` for `i` in `reservoir`. Note that these variables are always
nonnegative, as there is never incoming flow to a reservoir."
function variable_reservoir(wm::AbstractWaterModel, n::Int=wm.cnw)
    var(wm, n)[:qr] = JuMP.@variable(wm.model, [i in ids(wm, n, :reservoir)],
        base_name="qr[$(n)]", lower_bound=0.0,
        start=get_start(ref(wm, n, :reservoir), i, "qr_start", 0.0))
end

"Creates outgoing flow variables for all tanks in the network, i.e., `qt[i]`
for `i` in `tank`. Note that, unlike reservoirs, tanks can have inflow."
function variable_tank(wm::AbstractWaterModel, n::Int=wm.cnw)
    var(wm, n)[:qt] = JuMP.@variable(wm.model, [i in ids(wm, n, :tank)],
        base_name="qt[$(n)]",
        start=get_start(ref(wm, n, :node), i, "qt_start", 0.0))
end

"Creates bounded (by default) or unbounded (but still nonnegative) volume
variables for all tanks in the network, i.e., `V[i]` for `i` in `tank`."
function variable_volume(wm::AbstractWaterModel, n::Int=wm.cnw; bounded::Bool=true)
    bounded ? variable_volume_bounded(wm, n) : variable_volume_unbounded(wm, n)
end

"Creates bounded volume variables for all tanks in the network, i.e., `V[i]`
for `i` in `tank`. Bounds are related to minimum and maximum tank levels."
function variable_volume_bounded(wm::AbstractWaterModel, n::Int=wm.cnw)
    V_lb, V_ub = calc_tank_volume_bounds(wm, n)
    var(wm, n)[:V] = JuMP.@variable(wm.model, [i in ids(wm, n, :tank)],
        base_name="V[$(n)]", lower_bound=V_lb[i], upper_bound=V_ub[i],
        start=get_start(ref(wm, n, :tank), i, "V_start", V_ub[i]))
end

"Creates unbounded volume variables for all tanks in the network, i.e., `V[i]`
for `i` in `tank`. Note that since volume is inherently a nonnegative quantity,
however, these volume variables still contain lower bounds of zero."
function variable_volume_unbounded(wm::AbstractWaterModel, n::Int=wm.cnw)
    var(wm, n)[:V] = JuMP.@variable(wm.model, [i in ids(wm, n, :tank)],
        base_name="V[$(n)]", lower_bound=0.0,
        start=get_start(ref(wm, n, :tank), i, "V_start", 0.0))
end

### Link variables. ###
function variable_check_valve(wm::AbstractWaterModel, n::Int=wm.cnw)
    # Create variables for the statuses of check valves (one is open).
    var(wm, n)[:x_cv] = JuMP.@variable(wm.model,
        [i in ids(wm, n, :check_valve)], base_name="x_cv[$(n)]", binary=true,
        start=get_start(ref(wm, n, :link), i, "x_cv_start", 0))
end

function variable_flow_violation(wm::AbstractWaterModel, n::Int=wm.cnw)
    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    var(wm, n)[:Delta_e] = JuMP.@variable(wm.model, [a in ids(wm, n, :link)],
        base_name="Delta_e[$(n)]", lower_bound=0.0,
        start=get_start(ref(wm, n, :link), a, "Delta_e_start", 0.0))
end

function variable_fixed_speed_pump_operation(wm::AbstractWaterModel, n::Int=wm.cnw)
    var(wm, n)[:x_pump] = JuMP.@variable(wm.model, [a in ids(wm, n, :pump)],
        base_name="x_pump[$(n)]", binary=true,
        start=get_start(ref(wm, n, :pump), a, "x_pump_start", 1.0))
end

function variable_fixed_speed_pump_threshold(wm::AbstractWaterModel, n::Int=wm.cnw)
    var(wm, n)[:x_thrs_gt] = JuMP.@variable(wm.model, [a in ids(wm, n, :pump)],
        base_name="x_thrs_gt[$(n)]", binary=true,
        start=get_start(ref(wm, n, :pump), a, "x_thrs_gt_start", 0.0))
    var(wm, n)[:x_thrs_lt] = JuMP.@variable(wm.model, [a in ids(wm, n, :pump)],
        base_name="x_thrs_lt[$(n)]", binary=true,
        start=get_start(ref(wm, n, :pump), a, "x_thrs_lt_start", 0.0))
    var(wm, n)[:x_thrs_bt] = JuMP.@variable(wm.model, [a in ids(wm, n, :pump)],
        base_name="x_thrs_bt[$(n)]", binary=true,
        start=get_start(ref(wm, n, :pump), a, "x_thrs_bt_start", 1.0))
end

function variable_pump_operation(wm::AbstractWaterModel, n::Int=wm.cnw)
    variable_fixed_speed_pump_operation(wm, n)
end

function variable_pump_control(wm::AbstractWaterModel, n::Int=wm.cnw)
    variable_fixed_speed_pump_threshold(wm, n)
end

function variable_resistance(wm::AbstractWaterModel, n::Int=wm.cnw)
    link_ids = ids(wm, n, :link_ne)
    var(wm, n)[:x_res] = Dict{Int, Array{JuMP.VariableRef}}(a => [] for a in link_ids)

    for a in ids(wm, n, :link_ne)
        n_r = length(ref(wm, n, :resistance, a)) # Number of resistances.
        var(wm, n, :x_res)[a] = JuMP.@variable(wm.model, [r in 1:n_r], binary=true,
            lower_bound=0.0, upper_bound=1.0, base_name="x_res[$(n)][$(a)]",
            start=get_start(ref(wm, n, :link_ne), a, r, "x_res_start", 0.0))
    end
end
