# Define MILP-R implementations of water distribution models.

export CVXNLPWaterModel, StandardCVXNLPForm

import MathProgBase

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
    wm.var[:nw][n][:qp] = Dict{Int, Array{Variable, 1}}()
    wm.var[:nw][n][:qn] = Dict{Int, Array{Variable, 1}}()

    for (a, connection) in wm.ref[:nw][n][:connection]
        R_a = wm.ref[:nw][n][:resistance][a]

        # Initialize variables associated with flow from i to j.
        wm.var[:nw][n][:qp][a] = @variable(wm.model, [r in 1:length(R_a)],
                                           lowerbound = 0.0, start = 0.0,
                                           category = :Cont,
                                           basename = "qp_$(n)_$(a)")

        # Initialize variables associated with flow from j to i.
        wm.var[:nw][n][:qn][a] = @variable(wm.model, [r in 1:length(R_a)],
                                           lowerbound = 0.0, start = 0.0,
                                           category = :Cont,
                                           basename = "qn_$(n)_$(a)")
    end
end

#function get_head_solution(cvx::GenericWaterModel{T}, n::Int = cvx.cnw) where T <: StandardCVXNLPForm
#    junction_ids = collect(ids(cvx, n, :junctions))
#    reservoir_ids = collect(ids(cvx, n, :reservoirs))
#    connection_ids = collect(ids(cvx, n, :connection))
#    node_ids = [junction_ids; reservoir_ids]
#    node_mapping = Dict{Int, Int}(node_ids[i] => i for i in 1:length(node_ids))
#
#    num_reservoirs = length(reservoir_ids)
#    num_nodes = length(node_ids)
#    num_arcs = length(connection_ids)
#
#    # Create matrices for the left- and right-hand sides (for Ax = b).
#    A = zeros(Float64, num_arcs + num_reservoirs, num_nodes)
#    b = zeros(Float64, num_arcs + num_reservoirs, 1)
#
#    for (i, a) in enumerate(connection_ids)
#        connection = cvx.ref[:nw][n][:connection][a]
#
#        node_i = parse(Int, connection["node1"])
#        A[i, node_mapping[node_i]] = 1
#
#        node_j = parse(Int, connection["node2"])
#        A[i, node_mapping[node_j]] = -1
#
#        q_p = getvalue(cvx.var[:nw][n][:qp][a][1])
#        q_n = getvalue(cvx.var[:nw][n][:qn][a][1])
#
#        L = connection["length"]
#        resistance = cvx.ref[:nw][n][:resistance][a][1]
#
#        b[i] = L * resistance * (q_p - q_n) * abs(q_p - q_n)^(0.852)
#    end
#
#    for (i, reservoir_id) in enumerate(reservoir_ids)
#        A[num_arcs + i, node_mapping[reservoir_id]] = 1
#        b[num_arcs + i] = cvx.ref[:nw][n][:reservoirs][reservoir_id]["head"]
#    end
#
#    h = A \ b # Get the solution for head variables.
#    return Dict(i => h[node_mapping[i]] for i in node_ids)
#end
