########################################################################
# This file defines commonly-used variables for water systems models.
########################################################################

function get_start(set, item_key, value_key, default = 0.0)
    return get(get(set, item_key, Dict()), value_key, default)
end

function variable_undirected_flow(wm::GenericWaterModel{T}, n::Int=wm.cnw; bounded::Bool=true) where T <: AbstractWaterFormulation
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links)))

    # Initialize directed flow variables. The variables q⁺ correspond to flow
    # from i to j, and the variables q⁻ correspond to flow from j to i.
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
    wm.var[:nw][n][:qⁿᵉ] = Dict{Int, Array{JuMP.VariableRef, 1}}()

    # Initialize directed flow variables. The variables q⁺ correspond to flow
    # from i to j, and the variables q⁻ correspond to flow from j to i.
    if bounded
        lb, ub = calc_flow_rate_bounds(wm, n)

        for a in arcs
            num_resistances = length(wm.ref[:nw][n][:resistance][a])

            qⁿᵉ = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                 lower_bound=lb[a][r], upper_bound=ub[a][r],
                                 start=1.0e-6, base_name="qⁿᵉ[$(n)][$(a)]")

            wm.var[:nw][n][:qⁿᵉ][a] = qⁿᵉ
        end
    else
        for a in arcs
            num_resistances = length(wm.ref[:nw][n][:resistance][a])

            qⁿᵉ = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                 start = 1.0e-6, base_name = "qⁿᵉ[$(n)][$(a)]")

            wm.var[:nw][n][:qⁿᵉ][a] = qⁿᵉ
        end
    end
end

function variable_directed_flow(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852, bounded::Bool=true) where T <: AbstractWaterFormulation
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links)))

    # Initialize directed flow variables. The variables q⁺ correspond to flow
    # from i to j, and the variables q⁻ correspond to flow from j to i.
    if bounded
        ub⁻, ub⁺ = calc_directed_flow_upper_bounds(wm, alpha, n)

        q⁻ = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                            upper_bound=maximum(ub⁻[a]), base_name="q⁻[$(n)]",
                            start=get_start(wm.ref[:nw][n][:links], a,
                                            "q⁻_start", 1.0e-6))

        wm.var[:nw][n][:q⁻] = q⁻

        q⁺ = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                            upper_bound=maximum(ub⁺[a]), base_name="q⁺[$(n)]",
                            start=get_start(wm.ref[:nw][n][:links], a,
                                            "q⁺_start", 1.0e-6))

        wm.var[:nw][n][:q⁺] = q⁺
    else
        q⁻ = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                            base_name="q⁻[$(n)]",
                            start=get_start(wm.ref[:nw][n][:links], a,
                                            "q⁻_start", 1.0e-6))

        wm.var[:nw][n][:q⁻] = q⁻

        q⁺ = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                            base_name="q⁺[$(n)]",
                            start=get_start(wm.ref[:nw][n][:links], a,
                                            "q⁺_start", 1.0e-6))

        wm.var[:nw][n][:q⁺] = q⁺
    end
end

function variable_directed_flow_ne(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852, bounded::Bool=true) where T <: AbstractWaterFormulation
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links_ne)))
    resistances = wm.ref[:nw][n][:resistance]

    wm.var[:nw][n][:qⁿᵉ⁻] = Dict{Int, Array{JuMP.VariableRef, 1}}()
    wm.var[:nw][n][:qⁿᵉ⁺] = Dict{Int, Array{JuMP.VariableRef, 1}}()

    # Initialize directed flow variables. The variables q⁺ correspond to flow
    # from i to j, and the variables q⁻ correspond to flow from j to i.
    if bounded
        ub⁻, ub⁺ = calc_directed_flow_upper_bounds(wm, alpha, n)

        for a in arcs
            num_resistances = length(wm.ref[:nw][n][:resistance][a])

            qⁿᵉ⁻ = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                  lower_bound = 0.0, start = 0.0,
                                  upper_bound = max(0.0, ub⁻[a][r]),
                                  base_name = "qⁿᵉ⁻[$(n)][$(a)]")

            wm.var[:nw][n][:qⁿᵉ⁻][a] = qⁿᵉ⁻

            qⁿᵉ⁺ = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                  lower_bound = 0.0, start = 0.0,
                                  upper_bound = max(0.0, ub⁺[a][r]),
                                  base_name = "qⁿᵉ⁺[$(n)][$(a)]")

            wm.var[:nw][n][:qⁿᵉ⁺][a] = qⁿᵉ⁺
        end
    else
        for a in arcs
            num_resistances = length(wm.ref[:nw][n][:resistance][a])

            qⁿᵉ⁻ = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                  lower_bound = 0.0, start = 0.0,
                                  base_name = "qⁿᵉ⁻[$(n)][$(a)]")

            wm.var[:nw][n][:qⁿᵉ⁻][a] = qⁿᵉ⁻

            qⁿᵉ⁺ = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                  lower_bound = 0.0, start = 0.0,
                                  base_name = "qⁿᵉ⁺[$(n)][$(a)]")

            wm.var[:nw][n][:qⁿᵉ⁺][a] = qⁿᵉ⁺
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
                       upper_bound=ubs[i], start=ubs[i], base_name="h[$(n)]")
    wm.var[:nw][n][:h] = h
end

function variable_directed_head_difference(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links)))

    # Get the bounds for the head difference variables.
    lbs, ubs = calc_head_difference_bounds(wm, n)

    # Initialize variables associated with negative head differences.
    Δh⁻ = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                         upper_bound=abs(lbs[a]), start=abs(lbs[a]),
                         base_name="Δh⁻[$(n)]")
    wm.var[:nw][n][:Δh⁻] = Δh⁻

    # Initialize variables associated with positive head differences.
    Δh⁺ = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                         upper_bound=abs(ubs[a]), start=abs(ubs[a]),
                         base_name="Δh⁺[$(n)]")
    wm.var[:nw][n][:Δh⁺] = Δh⁺
end

function variable_flow_direction(wm::GenericWaterModel, n::Int=wm.cnw)
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links)))

    # Initialize variables associated with flow direction. If this variable is
    # equal to one, the flow direction is from i to j. If it is equal to zero,
    # the flow direction is from j to i.
    wm.var[:nw][n][:xᵈᶦʳ] = JuMP.@variable(wm.model, [a in arcs], start=1,
                                           binary=true, base_name="dir[$(n)]")
end

function variable_resistance_ne(wm::GenericWaterModel, n::Int=wm.cnw)
    # Get indices for all network arcs.
    arcs = sort(collect(ids(wm, n, :links)))

    # Initialize variables associated with flow direction. If this variable is
    # equal to one, the flow direction is from i to j. If it is equal to zero,
    # the flow direction is from j to i.
    wm.var[:nw][n][:xʳᵉˢ] = Dict{Int, Array{JuMP.VariableRef, 1}}()

    for (a, link) in wm.ref[:nw][n][:links]
        num_resistances = length(wm.ref[:nw][n][:resistance][a])

        wm.var[:nw][n][:xʳᵉˢ][a] = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                                 start=0, binary=true,
                                                 base_name="xʳᵉˢ[$(n)][$(a)]")

        JuMP.set_start_value(wm.var[:nw][n][:xʳᵉˢ][a][1], 1)
    end
end
