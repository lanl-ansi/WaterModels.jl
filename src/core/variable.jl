########################################################################
# This file defines commonly-used variables for water systems models.
########################################################################

function get_start(set, item_key, value_key, default = 0.0)
    return get(get(set, item_key, Dict()), value_key, default)
end

function variable_undirected_flow(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852, bounded::Bool=true) where T <: AbstractWaterFormulation
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n, :links))

    # Initialize directed flow variables. The variables q⁺ correspond to flow
    # from i to j, and the variables q⁻ correspond to flow from j to i.
    if bounded
        ub⁻, ub⁺ = calc_directed_flow_upper_bounds(wm, alpha, n)

        wm.var[:nw][n][:q] = JuMP.@variable(wm.model, [a in arcs],
                                            lower_bound = min(0.0, -maximum(ub⁻[a])),
                                            upper_bound = max(0.0, maximum(ub⁺[a])),
                                            start = 0.0,
                                            base_name = "q[$(n)]")
    else
        wm.var[:nw][n][:q] = JuMP.@variable(wm.model, [a in arcs],
                                            start = 0.0,
                                            base_name = "q[$(n)]")
    end
end

function variable_undirected_flow_ne(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852, bounded::Bool=true) where T <: AbstractWaterFormulation
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n, :links))
    resistances = wm.ref[:nw][n][:resistance]

    # Initialize directed flow variables. The variables q⁺ correspond to flow
    # from i to j, and the variables q⁻ correspond to flow from j to i.
    if bounded
        ub⁻, ub⁺ = calc_directed_flow_upper_bounds(wm, alpha, n)

        wm.var[:nw][n][:q_ne] = JuMP.@variable(wm.model, [a in arcs,
                                               r in 1:length(resistances[a])],
                                               lower_bound = min(0.0, -maximum(ub⁻[a])),
                                               upper_bound = max(0.0, maximum(ub⁺[a])),
                                               start = 0.0,
                                               base_name = "q_ne[$(n)]")
    else
        wm.var[:nw][n][:q_ne] = JuMP.@variable(wm.model, [a in arcs,
                                               r in 1:length(resistances[a])],
                                               start = 0.0,
                                               base_name = "q_ne[$(n)]")
    end

    println(wm.var[:nw][n][:q_ne])
end

function variable_directed_flow(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852, bounded::Bool=true) where T <: AbstractWaterFormulation
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n, :links))

    # Initialize directed flow variables. The variables q⁺ correspond to flow
    # from i to j, and the variables q⁻ correspond to flow from j to i.
    if bounded
        ub⁻, ub⁺ = calc_directed_flow_upper_bounds(wm, alpha, n)

        wm.var[:nw][n][:q⁻] = JuMP.@variable(wm.model, [a in arcs],
                                             lower_bound = 0.0,
                                             upper_bound = maximum(ub⁻[a]),
                                             start = 0.0,
                                             base_name = "q⁻[$(n)]")

        wm.var[:nw][n][:q⁺] = JuMP.@variable(wm.model, [a in arcs],
                                             lower_bound = 0.0,
                                             upper_bound = maximum(ub⁺[a]),
                                             start = 0.0,
                                             base_name = "q⁺[$(n)]")
    else
        wm.var[:nw][n][:q⁻] = JuMP.@variable(wm.model, [a in arcs],
                                             lower_bound = 0.0,
                                             start = 0.0,
                                             base_name = "q⁻[$(n)]")

        wm.var[:nw][n][:q⁺] = JuMP.@variable(wm.model, [a in arcs],
                                             lower_bound = 0.0,
                                             start = 0.0,
                                             base_name = "q⁺[$(n)]")
    end
end

function variable_directed_flow_ne(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852, bounded::Bool=true) where T <: AbstractWaterFormulation
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n, :links_ne))
    resistances = wm.ref[:nw][n][:resistance]

    # Initialize directed flow variables. The variables q⁺ correspond to flow
    # from i to j, and the variables q⁻ correspond to flow from j to i.
    if bounded
        ub⁻, ub⁺ = calc_directed_flow_upper_bounds(wm, alpha, n)

        wm.var[:nw][n][:q_ne⁻] = JuMP.@variable(wm.model, [a in arcs,
                                                r in 1:length(resistances[a])],
                                                lower_bound = 0.0,
                                                upper_bound = ub⁻[a][r],
                                                start = 0.0,
                                                base_name = "q_ne⁻[$(n)]")

        wm.var[:nw][n][:q_ne⁺] = JuMP.@variable(wm.model, [a in arcs,
                                                r in 1:length(resistances[a])],
                                                lower_bound = 0.0,
                                                upper_bound = ub⁺[a][r],
                                                start = 0.0,
                                                base_name = "q_ne⁺[$(n)]")
    else
        wm.var[:nw][n][:q_ne⁻] = JuMP.@variable(wm.model, [a in arcs,
                                                r in 1:length(resistances[a])],
                                                lower_bound = 0.0,
                                                start = 0.0,
                                                base_name = "q_ne⁻[$(n)]")

        wm.var[:nw][n][:q_ne⁺] = JuMP.@variable(wm.model, [a in arcs,
                                                r in 1:length(resistances[a])],
                                                lower_bound = 0.0,
                                                start = 0.0,
                                                base_name = "q_ne⁺[$(n)]")
    end
end

function variable_pressure_head(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Get indices for all network nodes.
    nodes = [collect(ids(wm, n, :junctions)); collect(ids(wm, n, :reservoirs))]

    # Get the bounds associated with heads at nodes.
    lbs, ubs = calc_head_bounds(wm, n)

    # Initialize variables associated with head.
    wm.var[:nw][n][:h] = JuMP.@variable(wm.model, [i in nodes],
                                        lower_bound = lbs[i],
                                        upper_bound = ubs[i],
                                        start = ubs[i],
                                        base_name = "h[$(n)]")
end

function variable_directed_head_difference(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractWaterFormulation
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n, :links))

    # Get the bounds for the head difference variables.
    lbs, ubs = calc_head_difference_bounds(wm, n)

    # Initialize variables associated with negative head differences.
    wm.var[:nw][n][:Δh⁻] = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                                          upper_bound=abs(lbs[a]), start=abs(lbs[a]),
                                          base_name="Δh⁻[$(n)]")

    # Initialize variables associated with positive head differences.
    wm.var[:nw][n][:Δh⁺] = JuMP.@variable(wm.model, [a in arcs], lower_bound=0.0,
                                          upper_bound=abs(ubs[a]), start=abs(ubs[a]),
                                          base_name="Δh⁺[$(n)]")
end

function variable_flow_direction(wm::GenericWaterModel, n::Int=wm.cnw)
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n, :links))

    # Initialize variables associated with flow direction. If this variable is
    # equal to one, the flow direction is from i to j. If it is equal to zero,
    # the flow direction is from j to i.
    wm.var[:nw][n][:dir] = JuMP.@variable(wm.model, [a in arcs], start=1,
                                          binary=true, base_name="dir[$(n)]")

    # Fix direction bounds if they can be fixed.
    dh_lbs, dh_ubs = calc_head_difference_bounds(wm, n)
    for (a, link) in wm.ref[:nw][n][:links]
        #if (dh_lbs[a] >= 0.0)
        #    JuMP.fix(wm.var[:nw][n][:dir][a], 1)
        #elseif (dh_ubs[a] <= 0.0)
        #    JuMP.fix(wm.var[:nw][n][:dir][a], 0)
        #end
    end
end

#function variable_undirected_flow(wm::GenericWaterModel, n::Int=wm.cnw)
#    # Get indices for all network arcs.
#    links = collect(ids(wm, n, :links))
#
#    wm.var[:nw][n][:q] = JuMP.@variable(wm.model, [a in links], start=0.0,
#                                        base_name="q[$(n)]")
#end

function variable_resistance_ne(wm::GenericWaterModel, n::Int=wm.cnw)
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n, :links))

    # Initialize variables associated with flow direction. If this variable is
    # equal to one, the flow direction is from i to j. If it is equal to zero,
    # the flow direction is from j to i.
    wm.var[:nw][n][:xr] = Dict{Int, Array{JuMP.VariableRef, 1}}()

    for (a, link) in wm.ref[:nw][n][:links]
        num_resistances = length(wm.ref[:nw][n][:resistance][a])

        wm.var[:nw][n][:xr][a] = JuMP.@variable(wm.model, [r in 1:num_resistances],
                                                start=0, binary=true,
                                                base_name="xr[$(n)][$(a)]")

        JuMP.set_start_value(wm.var[:nw][n][:xr][a][1], 1)
    end
end

##function variable_directed_flow(wm::GenericWaterModel, n::Int = wm.cnw) where T <: StandardMINLPForm
##    # Get indices for all network arcs.
##    arcs = collect(ids(wm, n, :link))
##
##    # Compute sets of resistances.
##    ub_n, ub_p = calc_directed_flow_upper_bounds(wm, n)
##
##    # Initialize directed flow variables. The variables qp correspond to flow
##    # from i to j, and the variables qn correspond to flow from j to i.
##    wm.var[:nw][n][:qp] = Dict{Int, Array{Variable, 1}}()
##    wm.var[:nw][n][:qn] = Dict{Int, Array{Variable, 1}}()
##
##    for (a, link) in wm.ref[:nw][n][:links]
##        R_a = wm.ref[:nw][n][:resistance][a]
##
##        # Initialize variables associated with flow from i to j.
##        wm.var[:nw][n][:qp][a] = @variable(wm.model, [r in 1:length(R_a)],
##                                           lower_bound = 0.0,
##                                           upper_bound = ub_p[a][r],
##                                           start = 0.0, category = :Cont,
##                                           base_name = "qp_$(n)_$(a)")
##
##        # Initialize flow for the variable with least resistance.
##        setvalue(wm.var[:nw][n][:qp][a][end], ub_p[a][end])
##
##        # Initialize variables associated with flow from j to i.
##        wm.var[:nw][n][:qn][a] = @variable(wm.model, [r in 1:length(R_a)],
##                                           lower_bound = 0.0,
##                                           upper_bound = ub_n[a][r],
##                                           start = 0.0, category = :Cont,
##                                           base_name = "qn_$(n)_$(a)")
##    end
##end
#
##function variable_flow(wm::GenericWaterModel, n::Int = wm.cnw)
##    lbs, ubs = calc_flow_bounds(wm.ref[:nw][n][:pipes])
##    wm.var[:nw][n][:q] = @variable(wm.model, [id in keys(wm.ref[:nw][n][:pipes])],
##                                   lower_bound = lbs[id], upper_bound = ubs[id],
##                                   base_name = "q_$(n)", start = 1.0e-4)
##end
##
##function variable_head_common(wm::GenericWaterModel, n::Int = wm.cnw)
##    lbs, ubs = calc_head_bounds(wm.ref[:nw][n][:junctions], wm.ref[:nw][n][:reservoirs])
##    junction_ids = [key for key in keys(wm.ref[:nw][n][:junctions])]
##    reservoir_ids = [key for key in keys(wm.ref[:nw][n][:reservoirs])]
##    ids = [junction_ids; reservoir_ids]
##
##    # Add the head variables to the model.
##    wm.var[:nw][n][:h] = @variable(wm.model, [i in ids], lower_bound = lbs[i],
##                                   upper_bound = ubs[i], base_name = "h_$(n)",
##                                   start = ubs[i])
##end
##
##"Variables associated with building pipes."
##function variable_pipe_ne_common(wm::GenericWaterModel, n::Int = wm.cnw)
##    # Set up required data to initialize junction variables.
##    pipe_ids = [key for key in keys(wm.ref[:nw][n][:ne_pipe])]
##    wm.var[:nw][n][:psi] = Dict{Int, Any}()
##    wm.var[:nw][n][:gamma] = Dict{Int, Any}()
##    wm.var[:nw][n][:gamma_sum] = Dict{Int, Any}()
##
##    for (pipe_id, pipe) in wm.ref[:nw][n][:ne_pipe]
##        diameters = [d["diameter"] for d in pipe["diameters"]]
##
##        # Create binary variables associated with whether or not a diameter is used.
##        wm.var[:nw][n][:psi][pipe_id] = @variable(wm.model, [d in diameters], category = :Bin,
##                                                  base_name = "psi_$(n)_$(pipe_id)", start = 0)
##        setvalue(wm.var[:nw][n][:psi][pipe_id][diameters[end]], 1)
##
##        # Create a variable that corresponds to the selection of lambda.
##        h_i = wm.var[:nw][n][:h][parse(Int, pipe["node1"])]
##        h_j = wm.var[:nw][n][:h][parse(Int, pipe["node2"])]
##
##        # Get additional data related to the variables.
##        hij_lb = getlower_bound(h_i) - getupper_bound(h_j)
##        hij_ub = getupper_bound(h_i) - getlower_bound(h_j)
##
##        if pipe["flow_direction"] == POSITIVE
##            hij_lb = max(0.0, hij_lb)
##        elseif pipe["flow_direction"] == NEGATIVE
##            hij_ub = min(0.0, hij_ub)
##        end
##
##        lbs = Dict(d => hij_lb / calc_friction_factor_hw_ne(pipe, d) for d in diameters)
##        ubs = Dict(d => hij_ub / calc_friction_factor_hw_ne(pipe, d) for d in diameters)
##
##        min_gamma, max_gamma = [minimum(collect(values(lbs))), maximum(collect(values(ubs)))]
##        wm.var[:nw][n][:gamma][pipe_id] = @variable(wm.model, [d in diameters],
##                                                    lower_bound = lbs[d], upper_bound = ubs[d],
##                                                    start = 0.5 * (ubs[d] + lbs[d]))
##        wm.var[:nw][n][:gamma_sum][pipe_id] = @variable(wm.model, lower_bound = min_gamma,
##                                                        upper_bound = max_gamma,
##                                                        start = 0.5 * (min_gamma + max_gamma))
##    end
##end
##
##function variable_objective_ne(wm::GenericWaterModel, n::Int = wm.cnw)
##    wm.var[:nw][n][:objective] = @variable(wm.model, base_name = "objective_$(n)",
##                                           lower_bound = 0.0, start = 0.0)
##end
