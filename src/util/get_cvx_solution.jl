export get_cvx_solution, check_solution_bounds

function get_cvx_solution(wm::GenericWaterModel,
                          resistance_indices::Dict{Int, Int},
                          solver::MathProgBase.AbstractMathProgSolver,
                          n::Int = wm.cnw)
    # Create a copy of the original network data.
    network = deepcopy(wm.data)

    for (a, connection) in wm.ref[:nw][n][:connection]
        # Remove unnecessary data from the network specification.
        delete!(network["pipes"][string(a)], "diameter")
        delete!(network["pipes"][string(a)], "diameters")
        delete!(network["pipes"][string(a)], "resistances")

        # Set the resistance to the selected resistance.
        resistances = wm.ref[:nw][n][:resistance][a]
        resistance_index = resistance_indices[a]
        selected_resistance = resistances[resistance_index]
        network["pipes"][string(a)]["resistance"] = selected_resistance
    end

    # Solve the CVXNLP model.
    cvx = build_generic_model(network, CVXNLPWaterModel, WaterModels.post_cvx_hw)
    setsolver(cvx.model, solver)
    status = JuMP.solve(cvx.model)

    # Get the solution.
    q = get_flow_solution(cvx)
    h = get_head_solution(cvx)

    # Return the solution.
    return q, h
end

function get_flow_solution(cvx::GenericWaterModel{T}, n::Int = cvx.cnw) where T <: StandardCVXNLPForm
    # Create a dictionary for the flow solution.
    connection_ids = collect(ids(cvx, n, :connection))
    q_sol = Dict{Int, Float64}(a => 0.0 for a in connection_ids)

    for (a, connection) in cvx.ref[:nw][n][:connection]
        qp_sol = getvalue(cvx.var[:nw][n][:qp][a][1])
        qn_sol = getvalue(cvx.var[:nw][n][:qn][a][1])
        q_sol[a] = qp_sol - qn_sol
    end

    return q_sol
end

function get_head_solution(cvx::GenericWaterModel{T}, n::Int = cvx.cnw) where T <: StandardCVXNLPForm
    junction_ids = collect(ids(cvx, n, :junctions))
    reservoir_ids = collect(ids(cvx, n, :reservoirs))
    connection_ids = collect(ids(cvx, n, :connection))
    node_ids = [junction_ids; reservoir_ids]
    node_mapping = Dict{Int, Int}(node_ids[i] => i for i in 1:length(node_ids))

    num_reservoirs = length(reservoir_ids)
    num_nodes = length(node_ids)
    num_arcs = length(connection_ids)

    # Create matrices for the left- and right-hand sides (for Ax = b).
    A = zeros(Float64, num_arcs + num_reservoirs, num_nodes)
    b = zeros(Float64, num_arcs + num_reservoirs, 1)

    for (row, a) in enumerate(connection_ids)
        connection = cvx.ref[:nw][n][:connection][a]

        node_i = parse(Int, connection["node1"])
        A[row, node_mapping[node_i]] = 1.0

        node_j = parse(Int, connection["node2"])
        A[row, node_mapping[node_j]] = -1.0

        q_p = getvalue(cvx.var[:nw][n][:qp][a][1])
        q_n = getvalue(cvx.var[:nw][n][:qn][a][1])

        L = connection["length"]
        resistance = cvx.ref[:nw][n][:resistance][a][1]
        b[row] = L * resistance * (q_p - q_n) * abs(q_p - q_n)^(0.852)
    end

    for (i, reservoir_id) in enumerate(reservoir_ids)
        A[num_arcs + i, node_mapping[reservoir_id]] = 1
        b[num_arcs + i] = cvx.ref[:nw][n][:reservoirs][reservoir_id]["head"]
    end

    h = A \ b # Get the solution for head variables.
    return Dict{Int, Float64}(i => h[node_mapping[i]] for i in node_ids)
end

function check_solution_bounds(wm::GenericWaterModel,
                               q::Dict{Int, Float64},
                               h::Dict{Int, Float64},
                               resistance_indices::Dict{Int, Int},
                               n::Int = wm.cnw)
    # Initialize dictionaries used to store arc infeasibility results.
    connection_ids = collect(ids(wm, n, :connection))
    q_sat_lb = Dict{Int, Bool}(a => true for a in connection_ids)
    q_sat_ub = Dict{Int, Bool}(a => true for a in connection_ids)

    # Compute bound satisfaction results for flow variables.
    for (a, connection) in wm.ref[:nw][n][:connection]
        # Get the selected resistance index for this arc.
        r_a = resistance_indices[a]

        if q[a] >= 0.0
            # Compute bound satisfaction for flow from i to j.
            q_sat_lb[a] = q[a] >= getlowerbound(wm.var[:nw][n][:qp][a][r_a])
            q_sat_ub[a] = q[a] <= getupperbound(wm.var[:nw][n][:qp][a][r_a])
        else
            # Compute bound satisfaction for flow from j to i.
            q_sat_lb[a] = -q[a] >= getlowerbound(wm.var[:nw][n][:qn][a][r_a])
            q_sat_ub[a] = -q[a] <= getupperbound(wm.var[:nw][n][:qn][a][r_a])
        end
    end

    # Initialize dictionaries used to store node infeasibility results.
    junction_ids = collect(ids(wm, n, :junctions))
    reservoir_ids = collect(ids(wm, n, :reservoirs))
    node_ids = [junction_ids; reservoir_ids]
    h_sat_lb = Dict{Int, Bool}(i => true for i in node_ids)
    h_sat_ub = Dict{Int, Bool}(i => true for i in node_ids)

    # Compute bound satisfaction results for head variables.
    for (i, junction) in wm.ref[:nw][n][:junctions]
        h_sat_lb[i] = h[i] >= getlowerbound(wm.var[:nw][n][:h][i])
        h_sat_ub[i] = h[i] <= getupperbound(wm.var[:nw][n][:h][i])
    end

    # Return dictionaries of variable bound satisfaction results.
    return q_sat_lb, q_sat_ub, h_sat_lb, h_sat_ub
end
