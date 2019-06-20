########################################################################
# This file defines commonly-used variables for water systems models.
########################################################################

function get_start(set, item_key, value_key, default)
    return get(get(set, item_key, Dict()), value_key, default)
end

function get_start(set, item_key, first_value_key, second_value_key, default)
    return get(get(get(set, item_key, Dict()), second_value_key, Dict()), first_value_key, default)
end

function variable_undirected_flow(wm::GenericWaterModel{T}, n::Int=wm.cnw; bounded::Bool=true) where T <: AbstractWaterFormulation
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links)))

    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    if bounded
        lb, ub = calc_flow_rate_bounds(wm, n)

        q = JuMP.@variable(wm.model, [a in arcs], lower_bound=minimum(lb[a]),
                           upper_bound=maximum(ub[a]), base_name="q[$(n)]",
                           start=get_start(wm.ref[:nw][n][:links], a,
                                           "q_start", 1.0e-6))

        wm.var[:nw][n][:q] = q
    else
        q = JuMP.@variable(wm.model, [a in arcs], base_name="q[$(n)]",
                           start=get_start(wm.ref[:nw][n][:links], a,
                                           "q_start", 1.0e-6))

        wm.var[:nw][n][:q] = q
    end
end

function variable_undirected_flow_ne(wm::GenericWaterModel{T}, n::Int=wm.cnw; bounded::Bool=true) where T <: AbstractWaterFormulation
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links)))
    wm.var[:nw][n][:q_ne] = Dict{Int, Array{JuMP.VariableRef, 1}}()

    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    if bounded
        lb, ub = calc_flow_rate_bounds(wm, n)

        for a in arcs
            num_resistances = length(wm.ref[:nw][n][:resistance][a])

            q_ne = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                 lower_bound=lb[a][r], upper_bound=ub[a][r],
                                 base_name="q_ne[$(n)][$(a)]",
                                 start=get_start(wm.ref[:nw][n][:links], a, r,
                                                 "q_ne_start", 1.0e-6))

            wm.var[:nw][n][:q_ne][a] = q_ne
        end
    else
        for a in arcs
            num_resistances = length(wm.ref[:nw][n][:resistance][a])

            q_ne = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                 base_name="q_ne[$(n)][$(a)]",
                                 start=get_start(wm.ref[:nw][n][:links], a, r,
                                                 "q_ne_start", 1.0e-6))

            wm.var[:nw][n][:q_ne][a] = q_ne
        end
    end
end

function variable_directed_flow(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852, bounded::Bool=true) where T <: AbstractWaterFormulation
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links)))

    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    if bounded
        ub_n, ub_p = calc_directed_flow_upper_bounds(wm, alpha, n)

        qn = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                            upper_bound=maximum(ub_n[a]), base_name="qn[$(n)]",
                            start=get_start(wm.ref[:nw][n][:links], a,
                                            "qn_start", 1.0e-6))

        wm.var[:nw][n][:qn] = qn

        qp = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                            upper_bound=maximum(ub_p[a]), base_name="qp[$(n)]",
                            start=get_start(wm.ref[:nw][n][:links], a,
                                            "qp_start", 1.0e-6))

        wm.var[:nw][n][:qp] = qp
    else
        qn = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                            base_name="qn[$(n)]",
                            start=get_start(wm.ref[:nw][n][:links], a,
                                            "qn_start", 1.0e-6))

        wm.var[:nw][n][:qn] = qn

        qp = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                            base_name="qp[$(n)]",
                            start=get_start(wm.ref[:nw][n][:links], a,
                                            "qp_start", 1.0e-6))

        wm.var[:nw][n][:qp] = qp
    end
end

function variable_directed_flow_ne(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852, bounded::Bool=true) where T <: AbstractWaterFormulation
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links_ne)))
    resistances = wm.ref[:nw][n][:resistance]

    wm.var[:nw][n][:qn_ne] = Dict{Int, Array{JuMP.VariableRef, 1}}()
    wm.var[:nw][n][:qp_ne] = Dict{Int, Array{JuMP.VariableRef, 1}}()

    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    if bounded
        ub_n, ub_p = calc_directed_flow_upper_bounds(wm, alpha, n)

        for a in arcs
            num_resistances = length(wm.ref[:nw][n][:resistance][a])

            qn_ne = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                  lower_bound = 0.0,
                                  upper_bound = max(0.0, ub_n[a][r]),
                                  start=get_start(wm.ref[:nw][n][:links_ne], a, r,
                                                  "qn_ne_start", 1.0e-6),
                                  base_name = "qn_ne[$(n)][$(a)]")

            wm.var[:nw][n][:qn_ne][a] = qn_ne

            qp_ne = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                  lower_bound = 0.0,
                                  upper_bound = max(0.0, ub_p[a][r]),
                                  start=get_start(wm.ref[:nw][n][:links_ne], a, r,
                                                  "qp_ne_start", 1.0e-6),
                                  base_name = "qp_ne[$(n)][$(a)]")

            wm.var[:nw][n][:qp_ne][a] = qp_ne
        end
    else
        for a in arcs
            num_resistances = length(wm.ref[:nw][n][:resistance][a])

            qn_ne = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                  lower_bound = 0.0,
                                  start=get_start(wm.ref[:nw][n][:links_ne], a, r,
                                                  "qn_ne_start", 1.0e-6),
                                  base_name = "qn_ne[$(n)][$(a)]")

            wm.var[:nw][n][:qn_ne][a] = qn_ne

            qp_ne = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                  lower_bound = 0.0,
                                  start=get_start(wm.ref[:nw][n][:links_ne], a, r,
                                                  "qp_ne_start", 1.0e-6),
                                  base_name = "qp_ne[$(n)][$(a)]")

            wm.var[:nw][n][:qp_ne][a] = qp_ne
        end
    end
end

function variable_pressure_head(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Get indices for all network nodes.
    junction_ids = sort(collect(ids(wm, n, :junctions)))

    # Get the bounds associated with heads at junctions.
    lbs, ubs = calc_head_bounds(wm, n)

    # Initialize variables associated with head.
    h = JuMP.@variable(wm.model, [i in junction_ids], lower_bound=lbs[i],
                       upper_bound=ubs[i], base_name="h[$(n)]",
                       start=get_start(wm.ref[:nw][n][:junctions], i,
                                       "h_start", ubs[i]))

    wm.var[:nw][n][:h] = h
end

function variable_directed_head_difference(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links)))

    # Get the bounds for the head difference variables.
    lbs, ubs = calc_head_difference_bounds(wm, n)

    # Initialize variables associated with negative head differences.
    dhn = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                         upper_bound=abs(lbs[a]),
                         start=get_start(wm.ref[:nw][n][:links], a,
                                         "dhn_start", 0.0),
                         base_name="dhn[$(n)]")

    wm.var[:nw][n][:dhn] = dhn

    # Initialize variables associated with positive head differences.
    dhp = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                         upper_bound=abs(ubs[a]),
                         start=get_start(wm.ref[:nw][n][:links], a,
                                         "dhp_start", 0.0),
                         base_name="dhp[$(n)]")

    wm.var[:nw][n][:dhp] = dhp
end

function variable_flow_direction(wm::GenericWaterModel, n::Int=wm.cnw)
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links)))

    # Initialize variables associated with flow direction. If this variable is
    # equal to one, the flow direction is from i to j. If it is equal to zero,
    # the flow direction is from j to i.
    x_dir = JuMP.@variable(wm.model, [a in arcs],
                          start=get_start(wm.ref[:nw][n][:links], a,
                                          "x_dir_start", 1.0),
                          binary=true, base_name="x_dir[$(n)]")

    wm.var[:nw][n][:x_dir] = x_dir
end

function variable_resistance_ne(wm::GenericWaterModel, n::Int=wm.cnw)
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links)))

    # Initialize variables associated with flow direction. If this variable is
    # equal to one, the flow direction is from i to j. If it is equal to zero,
    # the flow direction is from j to i.
    wm.var[:nw][n][:x_res] = Dict{Int, Array{JuMP.VariableRef, 1}}()

    for (a, link) in wm.ref[:nw][n][:links]
        num_resistances = length(wm.ref[:nw][n][:resistance][a])

        x_res = JuMP.@variable(wm.model, [r in 1:num_resistances], binary=true,
                              base_name="x_res[$(n)][$(a)]",
                              start=get_start(wm.ref[:nw][n][:links], a, r,
                                              "x_res_start", 0.0))

        wm.var[:nw][n][:x_res][a] = x_res
    end
end
