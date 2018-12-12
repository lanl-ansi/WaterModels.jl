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

function get_head_solution(wm::GenericWaterModel{T}, solver::MathProgBase.AbstractMathProgSolver, n::Int = wm.cnw) where T <: StandardCVXNLPForm
    # Initialie the model.
    model = Model(solver = solver)

    # Initialize variables associated with head.
    nodes = [collect(ids(wm, n, :junctions)); collect(ids(wm, n, :reservoirs))]
    h = @variable(model, [i in nodes], start = 0.0, category = :Cont)

    for (a, connection) in wm.ref[:nw][n][:connection]
        L = connection["length"]
        resistance = wm.ref[:nw][n][:resistance][a][1]

        qp = wm.var[:nw][n][:qp][a][1]
        qn = wm.var[:nw][n][:qn][a][1]
        q = getvalue(qp) - getvalue(qn)

        h_i = h[parse(Int, connection["node1"])]
        h_j = h[parse(Int, connection["node2"])]

        @constraint(model, sign(q) * h_i - h_j <= L * resistance * abs(q)^(1.852) + 1.0e-4)
        @constraint(model, sign(q) * h_i - h_j >= L * resistance * abs(q)^(1.852) - 1.0e-4)
    end

    for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
        @constraint(model, h[i] == reservoir["head"])
    end

    status = JuMP.solve(model)

    return getvalue(h)
end

#function solution_is_feasible(cvx::GenericWaterModel, wm::GenericWaterModel{T}, n::Int = wm.cnw) where T <: StandardCVXNLPForm
#    # Initialie the model.
#    model = Model(solver = solver)
#
#    # Initialize variables associated with head.
#    nodes = [collect(ids(wm, n, :junctions)); collect(ids(wm, n, :reservoirs))]
#    h = @variable(model, [i in nodes], start = 0.0, category = :Cont)
#
#    for (a, connection) in wm.ref[:nw][n][:connection]
#        L = connection["length"]
#        resistance = wm.ref[:nw][n][:resistance][a][1]
#
#        qp = cvx.var[:nw][n][:qp][a][1]
#        qn = cvx.var[:nw][n][:qn][a][1]
#
#        h_i = h[parse(Int, connection["node1"])]
#        h_j = h[parse(Int, connection["node2"])]
#
#        @constraint(model, h_i - h_j <= L * resistance * q * (q^2)^0.426 + 1.0e-4)
#        @constraint(model, h_i - h_j >= L * resistance * q * (q^2)^0.426 - 1.0e-4)
#    end
#
#    for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
#        @constraint(model, h[i] == reservoir["head"])
#    end
#
#    status = JuMP.solve(model)
#
#    return getvalue(h)
#end

#function constraint_select_resistance(wm::GenericWaterModel, a::Int, n::Int = wm.cnw)

#function variable_flow(wm::GenericWaterModel, n::Int = wm.cnw) where T <: AbstractCVXNLPForm
#    # Get all connections (e.g., pipes) in the network.
#    connections = wm.ref[:nw][n][:connection]
#
#    # Create the flow variables associated with positive directions (i to j).
#    wm.var[:nw][n][:qp] = @variable(wm.model, [id in keys(connections)],
#                                    category = :Cont, lowerbound = 0.0,
#                                    basename = "qp_$(n)", start = 1.0e-6)
#
#    # Create the flow variables associated with negative directions (j to i).
#    wm.var[:nw][n][:qn] = @variable(wm.model, [id in keys(connections)],
#                                    category = :Cont, lowerbound = 0.0,
#                                    basename = "qn_$(n)", start = 1.0e-6)
#end
#
#function constraint_flow_conservation(wm::GenericWaterModel, i::Int, n::Int = wm.cnw) where T <: AbstractCVXNLPForm
#    # Collect the required variables.
#    connections = wm.ref[:nw][n][:connection]
#
#    out_arcs = collect(keys(filter(is_out_node_function(i), connections)))
#    out_p_vars = Array{JuMP.Variable}([wm.var[:nw][n][:qp][a][1] for a in out_arcs])
#    out_n_vars = Array{JuMP.Variable}([wm.var[:nw][n][:qn][a][1] for a in out_arcs])
#
#    in_arcs = collect(keys(filter(is_in_node_function(i), connections)))
#    in_p_vars = Array{JuMP.Variable}([wm.var[:nw][n][:qp][a][1] for a in in_arcs])
#    in_n_vars = Array{JuMP.Variable}([wm.var[:nw][n][:qn][a][1] for a in in_arcs])
#
#    # Add the flow conservation constraints for junction nodes.
#    if !haskey(wm.con[:nw][n], :flow_conservation)
#        wm.con[:nw][n][:flow_conservation] = Dict{Int, ConstraintRef}()
#    end    
#
#    demand = wm.ref[:nw][n][:junctions][i]["demand"]
#
#    wm.con[:nw][n][:flow_conservation][i] = @constraint(wm.model, sum(in_p_vars) + sum(out_n_vars) -
#                                                        sum(out_p_vars) - sum(in_n_vars) == demand)
#end
#
function solution_is_feasible(wm::GenericWaterModel, n::Int = wm.cnw)
    num_junctions = length(wm.ref[:nw][n][:junctions])
    num_reservoirs = length(wm.ref[:nw][n][:reservoirs])
    num_nodes = num_junctions + num_reservoirs
    num_arcs = length(wm.ref[:nw][n][:connection])
    A = zeros(Float64, num_arcs + num_reservoirs, num_nodes)
    b = zeros(Float64, num_arcs + num_reservoirs, 1)

    for (a, connection) in wm.ref[:nw][n][:connection]
        A[a, parse(Int, connection["node1"])] = 1
        A[a, parse(Int, connection["node2"])] = -1
        q_p = getvalue(wm.var[:nw][n][:qp][a][1])
        q_n = getvalue(wm.var[:nw][n][:qn][a][1])
        resistance = connection["length"] * wm.ref[:nw][n][:resistance][a][1]
        b[a] = resistance * (q_p - q_n) * abs(q_p - q_n)^(0.852)
    end

    k = 1
    for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
        row = num_arcs + k
        A[row, i] = 1
        b[row] = reservoir["head"]
        k += 1
    end

    h = A \ b
    junctions = wm.ref[:nw][n][:junctions]
    reservoirs = wm.ref[:nw][n][:reservoirs]
    max_junc_elev = maximum([junc["elev"] for junc in values(junctions)])
    max_res_head = maximum([res["head"] for res in values(reservoirs)])
    max_head = max(max_junc_elev, max_res_head)

    for (i, junction) in wm.ref[:nw][n][:junctions]
        if h[i] < junction["elev"] || h[i] > max_head
            return false
        end
    end

    return true
end
