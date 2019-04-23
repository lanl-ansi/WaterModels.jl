# Define CVXNLP implementations of water distribution models.
export CVXNLPWaterModel, StandardCVXNLPForm

"AbstractCVXNLPForm is derived from AbstractWaterFormulation"
abstract type AbstractCVXNLPForm <: AbstractWaterFormulation end

"StandardCVXNLPForm is derived from AbstractCVXNLPForm"
abstract type StandardCVXNLPForm <: AbstractCVXNLPForm end

"The CVXNLP model relaxes constraints into the objective."
const CVXNLPWaterModel = GenericWaterModel{StandardCVXNLPForm}

"CVXNLP constructor."
CVXNLPWaterModel(data::Dict{String, Any}; kwargs...) = GenericWaterModel(data, StandardCVXNLPForm; kwargs...)

function variable_flow(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852) where T <: StandardCVXNLPForm
    variable_directed_flow(wm, n, alpha=alpha, bounded=false)
end

function variable_head(wm::GenericWaterModel{T}, n::Int=wm.cnw; alpha::Float64=1.852) where T <: StandardCVXNLPForm
end

function constraint_potential_loss(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: StandardCVXNLPForm
end

function constraint_flow_conservation(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: StandardCVXNLPForm
    constraint_directed_flow_conservation(wm, i, n)
end

#function variable_directed_flow(wm::GenericWaterModel{T}, n::Int=wm.cnw; bounded::Bool=true) where T <: StandardCVXNLPForm
#    # Get indices for all network arcs.
#    arcs = collect(ids(wm, n, :links))
#    ub⁻, ub⁺ = calc_directed_flow_upper_bounds(wm)
#
#    # Initialize directed flow variables. The variables q⁺ correspond to flow
#    # from i to j, and the variables q⁻ correspond to flow from j to i.
#    if bounded
#        wm.var[:nw][n][:q⁻] = JuMP.@variable(wm.model, [a in arcs],
#                                             lower_bound = 0.0,
#                                             upper_bound = maximum(ub⁻[a]),
#                                             start = 0.0,
#                                             base_name = "q⁻[$(n)]")
#
#        wm.var[:nw][n][:q⁺] = JuMP.@variable(wm.model, [a in arcs],
#                                             lower_bound = 0.0,
#                                             upper_bound = maximum(ub⁺[a]),
#                                             start = 0.0,
#                                             base_name = "q⁺[$(n)]")
#    else
#        #lbs, ubs = calc_flow_bounds(wm.ref[:nw][n][:pipes])
#        wm.var[:nw][n][:q⁻] = JuMP.@variable(wm.model, [a in arcs], lower_bound = 0.0,
#                                             start = 0.0, base_name = "q⁻[$(n)]")
#
#        wm.var[:nw][n][:q⁺] = JuMP.@variable(wm.model, [a in arcs], lower_bound = 0.0,
#                                         start = 0.0, base_name = "q⁺[$(n)]")
#    end
#end

#function constraint_undirected_flow_conservation(wm::GenericWaterModel{T}, i::Int, n::Int = wm.cnw) where T <: StandardCVXNLPForm
#    # Collect the required variables.
#    links = wm.ref[:nw][n][:links]
#
#    if !haskey(wm.con[:nw][n], :flow_conservation)
#        wm.con[:nw][n][:flow_conservation] = Dict{Int, JuMP.ConstraintRef}()
#    end
#
#    # Initialize the flow sum expression.
#    flow_sum = JuMP.@expression(wm.model, 0.0)
#
#    for (a, conn) in filter(a -> i == a.second["node2"], links)
#        flow_sum += wm.var[:nw][n][:q_p][a] - wm.var[:nw][n][:q_n][a]
#    end
#
#    for (a, conn) in filter(a -> i == a.second["node1"], links)
#        flow_sum += wm.var[:nw][n][:q_n][a] - wm.var[:nw][n][:q_p][a]
#    end
#
#    # Add the flow conservation constraint.
#    demand = wm.ref[:nw][n][:junctions][i]["demand"]
#    con = JuMP.@constraint(wm.model, flow_sum == demand)
#    wm.con[:nw][n][:flow_conservation][i] = con
#end
#
#function constraint_directed_flow_conservation(wm::GenericWaterModel{T}, i::Int, n::Int = wm.cnw) where T <: StandardCVXNLPForm
#    # Collect the required variables.
#    links = wm.ref[:nw][n][:links]
#
#    if !haskey(wm.con[:nw][n], :flow_conservation)
#        wm.con[:nw][n][:flow_conservation] = Dict{Int, JuMP.ConstraintRef}()
#    end
#
#    # Initialize the flow sum expression.
#    flow_sum = JuMP.@expression(wm.model, 0.0)
#
#    for (a, conn) in filter(a -> i == a.second["node2"], links)
#        flow_sum += wm.var[:nw][n][:q_p][a] - wm.var[:nw][n][:q_n][a]
#    end
#
#    for (a, conn) in filter(a -> i == a.second["node1"], links)
#        flow_sum += wm.var[:nw][n][:q_n][a] - wm.var[:nw][n][:q_p][a]
#    end
#
#    # Add the flow conservation constraint.
#    demand = wm.ref[:nw][n][:junctions][i]["demand"]
#    con = JuMP.@constraint(wm.model, flow_sum == demand)
#    wm.con[:nw][n][:flow_conservation][i] = con
#end

function objective_wf(wm::GenericWaterModel{T}, n::Int = wm.cnw) where T <: StandardCVXNLPForm
    links = wm.ref[:nw][n][:links]
    linear_expr = JuMP.@expression(wm.model, 0.0)

    for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
        for (a, link) in filter(a -> i == a.second["node1"], links)
            q⁻ = wm.var[:nw][n][:q⁻][a]
            q⁺ = wm.var[:nw][n][:q⁺][a]
            linear_expr -= reservoir["head"] * (q⁺ - q⁻)
        end

        for (a, link) in filter(a -> i == a.second["node2"], links)
            q⁻ = wm.var[:nw][n][:q⁻][a]
            q⁺ = wm.var[:nw][n][:q⁺][a]
            linear_expr -= reservoir["head"] * (q⁻ - q⁺)
        end
    end

    linear_term = JuMP.@variable(wm.model, base_name = "linear_objective_term")
    JuMP.@constraint(wm.model, linear_expr == linear_term)

    # Initialize the objective.
    objective_expr = JuMP.@NLexpression(wm.model, linear_term +
        sum(link["length"] * wm.ref[:nw][n][:resistance][a][1] *
        (if_alpha(wm.var[:nw][n][:q⁻][a]) + if_alpha(wm.var[:nw][n][:q⁺][a]))
        for (a, link) in wm.ref[:nw][n][:links]))

    return JuMP.@NLobjective(wm.model, MOI.MIN_SENSE, objective_expr)
end
