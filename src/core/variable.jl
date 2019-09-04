########################################################################
# This file defines commonly-used variables for water systems models.
########################################################################

function get_start(set, item_key, value_key, default)
    return get(get(set, item_key, Dict()), value_key, default)
end

function get_start(set, item_key, first_value_key, second_value_key, default)
    return get(get(get(set, item_key, Dict()), second_value_key, Dict()), first_value_key, default)
end

function variable_volume(wm::AbstractWaterModel, n::Int=wm.cnw)
    lb, ub = calc_tank_volume_bounds(wm, n)
    var(wm, n)[:V] = JuMP.@variable(wm.model, [i in ids(wm, n, :tanks)],
        base_name="V[$(n)]", lower_bound=lb[i], upper_bound=ub[i], 
        start=get_start(ref(wm, n, :nodes), i, "V_start", ub[i]))
end

function variable_check_valve(wm::AbstractWaterModel, n::Int=wm.cnw)
    var(wm, n)[:x_cv] = JuMP.@variable(wm.model, [i in ids(wm, n, :check_valves)],
        base_name="x_cv[$(n)]", binary=true,
        start=get_start(ref(wm, n, :links), i, "x_cv_start", 0.0))
end

function variable_undirected_flow(wm::AbstractWaterModel, n::Int=wm.cnw; bounded::Bool=true)
    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    if bounded
        lb, ub = calc_flow_rate_bounds(wm, n)
        var(wm, n)[:q] = JuMP.@variable(wm.model, [a in ids(wm, n, :links)],
            base_name="q[$(n)]", lower_bound=minimum(lb[a]),
            upper_bound=maximum(ub[a]),
            start=get_start(ref(wm, n, :links), a, "q_start", 1.0e-6))
    else
        var(wm, n)[:q] = JuMP.@variable(wm.model, [a in ids(wm, n, :links)],
            base_name="q[$(n)]",
            start=get_start(ref(wm, n, :links), a, "q_start", 1.0e-6))
    end
end

function variable_flow_violation(wm::AbstractWaterModel, n::Int=wm.cnw)
    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    var(wm, n)[:Delta_e] = JuMP.@variable(wm.model, [a in ids(wm, n, :links)],
        base_name="Delta_e[$(n)]", lower_bound=0.0,
        start=get_start(ref(wm, n, :links), a, "Delta_e_start", 0.0))
end

function variable_head_violation(wm::AbstractWaterModel, n::Int=wm.cnw)
    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    var(wm, n)[:Delta_v] = JuMP.@variable(wm.model, [i in ids(wm, n, :nodes)],
        base_name="Delta_v[$(n)]", lower_bound=0.0,
        start=get_start(ref(wm, n, :nodes), i, "Delta_v_start", 0.0))
end

function variable_undirected_flow(wm::AbstractDirectedFlowModel, n::Int=wm.cnw)
    var(wm, n)[:q] = JuMP.@expression(wm.model, [a in ids(wm, n, :links)],
        var(wm, n, :qp, a) - var(wm, n, :qn, a))
end

function variable_undirected_flow_ne(wm::AbstractWaterModel, n::Int=wm.cnw; bounded::Bool=true)
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links)))
    var(wm, n)[:q_ne] = Dict{Int, Array{JuMP.VariableRef, 1}}()

    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    if bounded
        lb, ub = calc_flow_rate_bounds(wm, n)

        for a in arcs
            num_resistances = length(ref(wm, n, :resistance, a))

            q_ne = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                 lower_bound=lb[a][r], upper_bound=ub[a][r],
                                 base_name="q_ne[$(n)][$(a)]",
                                 start=get_start(ref(wm, n, :links), a, r,
                                                 "q_ne_start", 1.0e-6))

            var(wm, n)[:q_ne][a] = q_ne
        end
    else
        for a in arcs
            num_resistances = length(ref(wm, n, :resistance, a))

            q_ne = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                 base_name="q_ne[$(n)][$(a)]",
                                 start=get_start(ref(wm, n, :links), a, r,
                                                 "q_ne_start", 1.0e-6))

            var(wm, n)[:q_ne][a] = q_ne
        end
    end
end

function variable_directed_flow(wm::AbstractWaterModel, n::Int=wm.cnw; bounded::Bool=true)
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links)))

    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    if bounded
        ub_n, ub_p = calc_directed_flow_upper_bounds(wm, ref(wm, n, :alpha), n)

        qp = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                            upper_bound=maximum(ub_p[a]), base_name="qp[$(n)]",
                            start=get_start(ref(wm, n, :links), a,
                                            "qp_start", 1.0e-6))

        var(wm, n)[:qp] = qp

        qn = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                            upper_bound=maximum(ub_n[a]), base_name="qn[$(n)]",
                            start=get_start(ref(wm, n, :links), a,
                                            "qn_start", 1.0e-6))

        var(wm, n)[:qn] = qn
    else
        qp = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                            base_name="qp[$(n)]",
                            start=get_start(ref(wm, n, :links), a,
                                            "qp_start", 1.0e-6))

        var(wm, n)[:qp] = qp

        qn = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                            base_name="qn[$(n)]",
                            start=get_start(ref(wm, n, :links), a,
                                            "qn_start", 1.0e-6))

        var(wm, n)[:qn] = qn
    end
end

function variable_directed_flow_ne(wm::AbstractWaterModel, n::Int=wm.cnw; bounded::Bool=true)
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links_ne)))
    resistances = ref(wm, n, :resistance)

    var(wm, n)[:qn_ne] = Dict{Int, Array{JuMP.VariableRef, 1}}()
    var(wm, n)[:qp_ne] = Dict{Int, Array{JuMP.VariableRef, 1}}()

    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    if bounded
        ub_n, ub_p = calc_directed_flow_upper_bounds(wm, ref(wm, n, :alpha), n)

        for a in arcs
            num_resistances = length(ref(wm, n, :resistance, a))

            qp_ne = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                  lower_bound = 0.0,
                                  upper_bound = max(0.0, ub_p[a][r]),
                                  start=get_start(ref(wm, n, :links_ne), a, r,
                                                  "qp_ne_start", 1.0e-6),
                                  base_name = "qp_ne[$(n)][$(a)]")

            var(wm, n)[:qp_ne][a] = qp_ne

            qn_ne = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                  lower_bound = 0.0,
                                  upper_bound = max(0.0, ub_n[a][r]),
                                  start=get_start(ref(wm, n, :links_ne), a, r,
                                                  "qn_ne_start", 1.0e-6),
                                  base_name = "qn_ne[$(n)][$(a)]")

            var(wm, n)[:qn_ne][a] = qn_ne
        end
    else
        for a in arcs
            num_resistances = length(ref(wm, n, :resistance, a))

            qp_ne = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                  lower_bound = 0.0,
                                  start=get_start(ref(wm, n, :links_ne), a, r,
                                                  "qp_ne_start", 1.0e-6),
                                  base_name = "qp_ne[$(n)][$(a)]")

            var(wm, n)[:qp_ne][a] = qp_ne

            qn_ne = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                  lower_bound = 0.0,
                                  start=get_start(ref(wm, n, :links_ne), a, r,
                                                  "qn_ne_start", 1.0e-6),
                                  base_name = "qn_ne[$(n)][$(a)]")

            var(wm, n)[:qn_ne][a] = qn_ne
        end
    end
end

function variable_hydraulic_head(wm::AbstractWaterModel, n::Int=wm.cnw; bounded::Bool=true)
    if bounded
        lb, ub = calc_head_bounds(wm, n)
        var(wm, n)[:h] = JuMP.@variable(wm.model, [i in ids(wm, n, :nodes)],
            base_name="h[$(n)]", lower_bound=lb[i], upper_bound=ub[i],
            start=get_start(ref(wm, n, :nodes), i, "h_start", ub[i]))
    else
        var(wm, n)[:h] = JuMP.@variable(wm.model, [i in ids(wm, n, :nodes)],
            base_name="h[$(n)]",
            start=get_start(ref(wm, n, :nodes), i, "h_start", 0.0))
    end
end

function variable_directed_head_difference(wm::AbstractWaterModel, n::Int=wm.cnw)
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links)))

    # Get the bounds for the head difference variables.
    lbs, ubs = calc_head_difference_bounds(wm, n)

    # Initialize variables associated with positive head differences.
    dhp = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                         upper_bound=abs(ubs[a]),
                         start=get_start(ref(wm, n, :links), a,
                                         "dhp_start", 0.0),
                         base_name="dhp[$(n)]")

    var(wm, n)[:dhp] = dhp

    # Initialize variables associated with negative head differences.
    dhn = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                         upper_bound=abs(lbs[a]),
                         start=get_start(ref(wm, n, :links), a,
                                         "dhn_start", 0.0),
                         base_name="dhn[$(n)]")

    var(wm, n)[:dhn] = dhn
end

function variable_flow_direction(wm::AbstractWaterModel, n::Int=wm.cnw)
    # Initialize variables associated with flow direction. If this variable is
    # equal to one, the flow direction is from i to j. If it is equal to zero,
    # the flow direction is from j to i.
    var(wm, n)[:x_dir] = JuMP.@variable(wm.model, [a in ids(wm, n, :links)],
        base_name="x_dir[$(n)]", binary=true,
        start=get_start(ref(wm, n, :links), a, "x_dir_start", 1.0))
end

function variable_fixed_speed_pump_operation(wm::AbstractWaterModel, n::Int=wm.cnw)
    var(wm, n)[:x_pump] = JuMP.@variable(wm.model, [a in ids(wm, n, :pumps)],
        base_name="x_pump[$(n)]", binary=true,
        start=get_start(ref(wm, n, :pumps), a, "x_pump_start", 1.0))
end

function variable_fixed_speed_pump_threshold(wm::AbstractWaterModel, n::Int=wm.cnw)
    var(wm, n)[:x_thrs_gt] = JuMP.@variable(wm.model, [a in ids(wm, n, :pumps)],
        base_name="x_thrs_gt[$(n)]", binary=true,
        start=get_start(ref(wm, n, :pumps), a, "x_thrs_gt_start", 0.0))
    var(wm, n)[:x_thrs_lt] = JuMP.@variable(wm.model, [a in ids(wm, n, :pumps)],
        base_name="x_thrs_lt[$(n)]", binary=true,
        start=get_start(ref(wm, n, :pumps), a, "x_thrs_lt_start", 0.0))
    var(wm, n)[:x_thrs_bt] = JuMP.@variable(wm.model, [a in ids(wm, n, :pumps)],
        base_name="x_thrs_bt[$(n)]", binary=true,
        start=get_start(ref(wm, n, :pumps), a, "x_thrs_bt_start", 1.0))
end

function variable_head_gain(wm::AbstractWaterModel, n::Int=wm.cnw)
    var(wm, n)[:g] = JuMP.@variable(wm.model, [a in ids(wm, n, :pumps)],
        base_name="g[$(n)]", lower_bound=0.0,
        start=get_start(ref(wm, n, :links), a, "g", 0.0))
end

function variable_resistance_ne(wm::AbstractWaterModel, n::Int=wm.cnw)
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links)))

    # Initialize variables associated with flow direction. If this variable is
    # equal to one, the flow direction is from i to j. If it is equal to zero,
    # the flow direction is from j to i.
    var(wm, n)[:x_res] = Dict{Int, Array{JuMP.VariableRef, 1}}()

    for (a, link) in ref(wm, n, :links)
        num_resistances = length(ref(wm, n, :resistance, a))

        x_res = JuMP.@variable(wm.model, [r in 1:num_resistances], binary=true,
                              base_name="x_res[$(n)][$(a)]",
                              start=get_start(ref(wm, n, :links), a, r,
                                              "x_res_start", 0.0))

        var(wm, n)[:x_res][a] = x_res
    end
end

function variable_reservoir(wm::AbstractWaterModel, n::Int=wm.cnw)
    # Initialize variables associated with reservoir flow.
    var(wm, n)[:q_r] = JuMP.@variable(wm.model, [i in ids(wm, n, :reservoirs)],
        base_name="q_r[$(n)]",
        start=get_start(ref(wm, n, :reservoirs), i, "q_r_start", 0.0))
end

function variable_tank(wm::AbstractWaterModel, n::Int=wm.cnw)
    var(wm, n)[:q_t] = JuMP.@variable(wm.model, [i in ids(wm, n, :tanks)],
        base_name="q_t[$(n)]",
        start=get_start(ref(wm, n, :nodes), i, "q_t_start", 0.0))
end
