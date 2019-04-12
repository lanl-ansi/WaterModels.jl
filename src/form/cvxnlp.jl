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

function variable_directed_flow(wm::GenericWaterModel{T}, n::Int = wm.cnw) where T <: StandardCVXNLPForm
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n, :connection))

    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    wm.var[:nw][n][:q_n] = JuMP.@variable(wm.model, [a in arcs], lower_bound = 0.0,
                                          start = 0.0, base_name = "q_n[$(n)]")

    wm.var[:nw][n][:q_p] = JuMP.@variable(wm.model, [a in arcs], lower_bound = 0.0,
                                          start = 0.0, base_name = "q_p[$(n)]")
end

function constraint_directed_flow_conservation(wm::GenericWaterModel{T}, i::Int, n::Int = wm.cnw) where T <: StandardCVXNLPForm
    # Collect the required variables.
    connections = wm.ref[:nw][n][:connection]

    if !haskey(wm.con[:nw][n], :flow_conservation)
        wm.con[:nw][n][:flow_conservation] = Dict{Int, JuMP.ConstraintRef}()
    end

    # Initialize the flow sum expression.
    flow_sum = JuMP.@expression(wm.model, 0.0)

    for (a, conn) in filter(a -> i == parse(Int, a.second["node2"]), connections)
        flow_sum += wm.var[:nw][n][:q_p][a] - wm.var[:nw][n][:q_n][a]
    end

    for (a, conn) in filter(a -> i == parse(Int, a.second["node1"]), connections)
        flow_sum += wm.var[:nw][n][:q_n][a] - wm.var[:nw][n][:q_p][a]
    end

    # Add the flow conservation constraint.
    demand = wm.ref[:nw][n][:junctions][i]["demand"]
    con = JuMP.@constraint(wm.model, flow_sum == demand)
    wm.con[:nw][n][:flow_conservation][i] = con
end

function objective_wf(wm::GenericWaterModel{T}, n::Int = wm.cnw) where T <: StandardCVXNLPForm
    connections = wm.ref[:nw][n][:connection]
    linear_expr = JuMP.@expression(wm.model, 0.0)

    for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
        for (a, connection) in filter(a -> i == parse(Int, a.second["node1"]), connections)
            q_n = wm.var[:nw][n][:q_n][a]
            q_p = wm.var[:nw][n][:q_p][a]
            linear_expr -= reservoir["head"] * (q_p - q_n)
        end

        for (a, connection) in filter(a -> i == parse(Int, a.second["node2"]), connections)
            q_n = wm.var[:nw][n][:q_n][a]
            q_p = wm.var[:nw][n][:q_p][a]
            linear_expr -= reservoir["head"] * (q_n - q_p)
        end
    end

    linear_term = JuMP.@variable(wm.model, base_name = "linear_objective_term")
    JuMP.@constraint(wm.model, linear_expr == linear_term)

    # Initialize the objective.
    objective_expr = JuMP.@NLexpression(wm.model, linear_term +
        sum(connection["length"] * wm.ref[:nw][n][:resistance][a][1] *
        (if_alpha(wm.var[:nw][n][:q_n][a]) + if_alpha(wm.var[:nw][n][:q_p][a]))
        for (a, connection) in wm.ref[:nw][n][:connection]))

    return JuMP.@NLobjective(wm.model, MOI.MIN_SENSE, objective_expr)
end
