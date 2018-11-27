# Define MILP-R implementations of water distribution models.

export CVXNLPWaterModel, StandardCVXNLPForm

"AbstractCVXNLPForm is derived from AbstractWaterFormulation"
abstract type AbstractCVXNLPForm <: AbstractWaterFormulation end

"StandardCVXNLPForm is derived from AbstractCVXNLPForm"
abstract type StandardCVXNLPForm <: AbstractCVXNLPForm end

"The CVXNLP model relaxes constraints into the objective."
const CVXNLPWaterModel = GenericWaterModel{StandardCVXNLPForm}

"CVXNLP constructor."
CVXNLPWaterModel(data::Dict{String, Any}; kwargs...) = GenericWaterModel(data, StandardCVXNLPForm; kwargs...)

function variable_flow_cvxnlp(wm::GenericWaterModel, n::Int = wm.cnw) where T <: AbstractCVXNLPForm
    # Get all connections (e.g., pipes) in the network.
    connections = wm.ref[:nw][n][:connection]

    # Create the flow variables associated with positive directions (i to j).
    wm.var[:nw][n][:qp] = @variable(wm.model, [id in keys(connections)],
                                    category = :Cont, lowerbound = 0.0,
                                    basename = "qp_$(n)", start = 1.0e-6)

    # Create the flow variables associated with negative directions (j to i).
    wm.var[:nw][n][:qn] = @variable(wm.model, [id in keys(connections)],
                                    category = :Cont, lowerbound = 0.0,
                                    basename = "qn_$(n)", start = 1.0e-6)
end

function constraint_flow_conservation_cvx(wm::GenericWaterModel, i::Int, n::Int = wm.cnw) where T <: AbstractCVXNLPForm
    # Collect the required variables.
    connections = wm.ref[:nw][n][:connection]

    out_arcs = collect(keys(filter(is_out_node_function(i), connections)))
    out_p_vars = Array{JuMP.Variable}([wm.var[:nw][n][:qp][a] for a in out_arcs])
    out_n_vars = Array{JuMP.Variable}([wm.var[:nw][n][:qn][a] for a in out_arcs])

    in_arcs = collect(keys(filter(is_in_node_function(i), connections)))
    in_p_vars = Array{JuMP.Variable}([wm.var[:nw][n][:qp][a] for a in in_arcs])
    in_n_vars = Array{JuMP.Variable}([wm.var[:nw][n][:qn][a] for a in in_arcs])

    # Add the flow conservation constraints for junction nodes.
    if !haskey(wm.con[:nw][n], :flow_conservation)
        wm.con[:nw][n][:flow_conservation] = Dict{Int, ConstraintRef}()
    end    

    demand = wm.ref[:nw][n][:junctions][i]["demand"]

    wm.con[:nw][n][:flow_conservation][i] = @constraint(wm.model, sum(in_p_vars) + sum(out_n_vars) -
                                                        sum(out_p_vars) - sum(in_n_vars) == demand)
end

function solution_is_feasible(wm::GenericWaterModel{T}, n::Int = wm.cnw) where T <: AbstractCVXNLPForm
    num_junctions = length(wm.ref[:nw][n][:junctions])
    num_reservoirs = length(wm.ref[:nw][n][:reservoirs])
    num_nodes = num_junctions + num_reservoirs
    num_arcs = length(wm.ref[:nw][n][:connection])
    A = zeros(Float64, num_arcs + num_reservoirs, num_nodes)
    b = zeros(Float64, num_arcs + num_reservoirs, 1)

    for (a, connection) in wm.ref[:nw][n][:connection]
        A[a, parse(Int, connection["node1"])] = 1
        A[a, parse(Int, connection["node2"])] = -1
        q_p = getvalue(wm.var[:nw][n][:qp][a])
        q_n = getvalue(wm.var[:nw][n][:qn][a])
        r = calc_resistance_per_length_hw(connection)
        q = (q_p - q_n)
        b[a] = connection["length"] * r * q * abs(q)^(0.852)
    end

    k = 1
    for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
        row = num_arcs + k
        A[row, i] = 1
        b[row] = reservoir["head"]
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
