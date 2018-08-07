########################################################################
# This file defines commonly-used constraints for water systems models.
########################################################################

function construct_hw_separators(q::JuMP.Variable, lambda::Float64, n::Int = 5)
    q_points = linspace(getlowerbound(q), getupperbound(q), n)
    f_evals = [lambda * (x^2)^0.926 for x in q_points]
    df_evals = [lambda * 1.852*x / (x^2)^0.074 for x in q_points]
    df_evals[isnan.(df_evals)] = 0.0
    return [f_evals[i] + (q - q_points[i]) * df_evals[i] for i in 1:n]
end

function construct_dw_separators(q::JuMP.Variable, lambda::Float64, n::Int = 5)
    q_points = linspace(getlowerbound(q), getupperbound(q), n)
    f_evals = [lambda * x^2 for x in q_points]
    df_evals = [2 * lambda * x for x in q_points]
    return [f_evals[i] + (q - q_points[i]) * df_evals[i] for i in 1:n]
end

function constraint_flow_conservation_unknown_directions{T}(wm::GenericWaterModel{T}, i, n::Int = wm.cnw)
    # Add the flow conservation constraints for junction nodes.
    out_arcs = collect(keys(filter((id, pipe) -> pipe["node1"] == i, wm.ref[:nw][n][:pipes])))
    out_vars = Array{JuMP.Variable}([wm.var[:nw][n][:q][a] for a in out_arcs])
    in_arcs = collect(keys(filter((id, pipe) -> pipe["node2"] == i, wm.ref[:nw][n][:pipes])))
    in_vars = Array{JuMP.Variable}([wm.var[:nw][n][:q][a] for a in in_arcs])

    if haskey(wm.ref[:nw][n][:junctions], i)
        demand = wm.ref[:nw][n][:junctions][i]["demand"]
        @constraint(wm.model, sum(in_vars) - sum(out_vars) == demand)
    elseif haskey(wm.ref[:nw][n][:reservoirs], i)
        junctions = values(wm.ref[:nw][n][:junctions])
        sum_demand = sum(junction["demand"] for junction in junctions)
        @constraint(wm.model, sum(out_vars) - sum(in_vars) >= 0.0)
        @constraint(wm.model, sum(out_vars) - sum(in_vars) <= sum_demand)
    end
end

function constraint_flow_conservation_known_directions{T}(wm::GenericWaterModel{T}, i, n::Int = wm.cnw)
    # Add the flow conservation constraints for junction nodes.
    out_arcs = collect(keys(filter((id, pipe) -> pipe["node1"] == i, wm.ref[:nw][n][:pipes])))
    out_vars = Array{JuMP.Variable}([wm.var[:nw][n][:q][a] for a in out_arcs])
    out_dirs = [Int(wm.data["pipes"][a]["flow_direction"]) for a in out_arcs]
    flow_out = sum(dot(out_vars, out_dirs))

    in_arcs = collect(keys(filter((id, pipe) -> pipe["node2"] == i, wm.ref[:nw][n][:pipes])))
    in_vars = Array{JuMP.Variable}([wm.var[:nw][n][:q][a] for a in in_arcs])
    in_dirs = [Int(wm.data["pipes"][a]["flow_direction"]) for a in in_arcs]
    flow_in = sum(dot(in_vars, in_dirs))

    if haskey(wm.ref[:nw][n][:junctions], i)
        demand = wm.ref[:nw][n][:junctions][i]["demand"]
        @constraint(wm.model, flow_in - flow_out == demand)
    elseif haskey(wm.ref[:nw][n][:reservoirs], i)
        junctions = values(wm.ref[:nw][n][:junctions])
        sum_demand = sum(junction["demand"] for junction in junctions)
        @constraint(wm.model, flow_out - flow_in >= 0.0)
        @constraint(wm.model, flow_out - flow_in <= sum_demand)
    end
end

function constraint_potential_flow_coupling_relaxed_hw{T}(wm::GenericWaterModel{T}, a,
                                                          num_separators::Int = 5,
                                                          n::Int = wm.cnw)
    gamma = wm.var[:nw][n][:gamma][a]
    q = wm.var[:nw][n][:q][a]
    lambda = calc_friction_factor_hw(wm.ref[:nw][n][:pipes][a])

    # Use the piecewise linear outer approximation.
    for cut in construct_hw_separators(q, lambda, num_separators)
        @constraint(wm.model, gamma >= cut)
    end
end

function constraint_potential_flow_coupling_exact_hw{T}(wm::GenericWaterModel{T}, a,
                                                        n::Int = wm.cnw)
    gamma = wm.var[:nw][n][:gamma][a]
    q = wm.var[:nw][n][:q][a]
    lambda = calc_friction_factor_hw(wm.ref[:nw][n][:pipes][a])

    setlowerbound(q, 0.0) # In this version of the problem, this variable is nonnegative.
    @NLconstraint(wm.model, gamma == lambda * (q + 1e-7)^1.852)
end

function constraint_potential_flow_coupling_relaxed_dw{T}(wm::GenericWaterModel{T}, a,
                                                          num_separators::Int = 5,
                                                          n::Int = wm.cnw)
    gamma = wm.var[:nw][n][:gamma][a]
    q = wm.var[:nw][n][:q][a]
    viscosity = wm.ref[:nw][n][:options]["viscosity"]
    lambda = calc_friction_factor_dw(wm.ref[:nw][n][:pipes][a], viscosity)

    # Add the constraint.
    @constraint(wm.model, gamma >= lambda * q^2)

    # # Use the piecewise linear outer approximation.
    # for cut in construct_dw_separators(q, lambda, num_separators)
    #     @constraint(wm.model, gamma >= cut)
    # end
end

function constraint_potential_flow_coupling_exact_dw{T}(wm::GenericWaterModel{T}, a,
                                                        n::Int = wm.cnw)
    gamma = wm.var[:nw][n][:gamma][a]
    q = wm.var[:nw][n][:q][a]
    viscosity = wm.ref[:nw][n][:options]["viscosity"]
    lambda = calc_friction_factor_dw(wm.ref[:nw][n][:pipes][a], viscosity)

    # Add the full constraint.
    @NLconstraint(wm.model, gamma == lambda * q^2)
end

function constraint_define_gamma_unknown_direction{T}(wm::GenericWaterModel{T},
                                                      a, n::Int = wm.cnw)
    # Get source and target nodes corresponding to the arc.
    i = wm.ref[:nw][n][:pipes][a]["node1"]
    j = wm.ref[:nw][n][:pipes][a]["node2"]

    # Collect variables needed for the constraint.
    gamma = wm.var[:nw][n][:gamma][a]
    y_p = wm.var[:nw][n][:yp][a]
    y_n = wm.var[:nw][n][:yn][a]
    h_i = wm.var[:nw][n][:h][i]
    h_j = wm.var[:nw][n][:h][j]

    # Get additional data related to the variables.
    h_i_lb = getlowerbound(h_i)
    h_i_ub = getupperbound(h_i)
    h_j_lb = getlowerbound(h_j)
    h_j_ub = getupperbound(h_j)

    # Add the required constraints.
    @constraint(wm.model, h_j - h_i + (h_i_lb - h_j_ub) * (y_p - y_n + 1) <= gamma)
    @constraint(wm.model, h_i - h_j + (h_i_ub - h_j_lb) * (y_p - y_n - 1) <= gamma)
    @constraint(wm.model, h_j - h_i + (h_i_ub - h_j_lb) * (y_p - y_n + 1) >= gamma)
    @constraint(wm.model, h_i - h_j + (h_i_lb - h_j_ub) * (y_p - y_n - 1) >= gamma)
end

function constraint_define_gamma_negative_direction{T}(wm::GenericWaterModel{T},
                                                       a, n::Int = wm.cnw)
    # Get source and target nodes corresponding to the arc.
    i = wm.ref[:nw][n][:pipes][a]["node1"]
    j = wm.ref[:nw][n][:pipes][a]["node2"]

    # Add the required constraint.
    gamma = wm.var[:nw][n][:gamma][a]
    h_i = wm.var[:nw][n][:h][i]
    h_j = wm.var[:nw][n][:h][j]
    @constraint(wm.model, h_j - h_i == gamma)
end

function constraint_define_gamma_positive_direction{T}(wm::GenericWaterModel{T},
                                                       a, n::Int = wm.cnw)
    # Get source and target nodes corresponding to the arc.
    i = wm.ref[:nw][n][:pipes][a]["node1"]
    j = wm.ref[:nw][n][:pipes][a]["node2"]

    # Add the required constraint.
    gamma = wm.var[:nw][n][:gamma][a]
    h_i = wm.var[:nw][n][:h][i]
    h_j = wm.var[:nw][n][:h][j]
    @constraint(wm.model, h_i - h_j == gamma)
end

function constraint_bidirectional_flow{T}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Get source and target nodes corresponding to the arc.
    i = wm.ref[:nw][n][:pipes][a]["node1"]
    j = wm.ref[:nw][n][:pipes][a]["node2"]

    # Collect variables needed for the constraint.
    q = wm.var[:nw][n][:q][a]
    y_p = wm.var[:nw][n][:yp][a]
    y_n = wm.var[:nw][n][:yn][a]
    h_i = wm.var[:nw][n][:h][i]
    h_j = wm.var[:nw][n][:h][j]

    # Get additional data related to the variables.
    h_i_lb = getlowerbound(h_i)
    h_i_ub = getupperbound(h_i)
    h_j_lb = getlowerbound(h_j)
    h_j_ub = getupperbound(h_j)

    # Get the sum of all junction demands.
    junctions = values(wm.ref[:nw][n][:junctions])
    sum_demand = sum(junction["demand"] for junction in junctions)

    # Add the first set of constraints.
    @constraint(wm.model, (y_p - 1) * sum_demand <= q)
    @constraint(wm.model, (1 - y_n) * sum_demand >= q)

    # Add the second set of constraints.
    @constraint(wm.model, (1 - y_p) * (h_i_lb - h_j_ub) <= h_i - h_j)
    @constraint(wm.model, (1 - y_n) * (h_i_ub - h_j_lb) >= h_i - h_j)

    # Add the third set of constraints.
    @constraint(wm.model, y_p + y_n == 1)
end

function constraint_negative_flow{T}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Get source and target nodes corresponding to the arc.
    i = wm.ref[:nw][n][:pipes][a]["node1"]
    j = wm.ref[:nw][n][:pipes][a]["node2"]

    # Collect variables and constants needed for the constraint.
    q = wm.var[:nw][n][:q][a]
    junctions = values(wm.ref[:nw][n][:junctions])
    sum_demand = sum(junction["demand"] for junction in junctions)

    # Add the necessary constraints.
    @constraint(wm.model, q <= 0.0)
    @constraint(wm.model, q >= -sum_demand)
end

function constraint_positive_flow{T}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Get source and target nodes corresponding to the arc.
    i = wm.ref[:nw][n][:pipes][a]["node1"]
    j = wm.ref[:nw][n][:pipes][a]["node2"]

    # Collect variables and constants needed for the constraint.
    q = wm.var[:nw][n][:q][a]
    junctions = values(wm.ref[:nw][n][:junctions])
    sum_demand = sum(junction["demand"] for junction in junctions)

    # Add the necessary constraints.
    @constraint(wm.model, q >= 0.0)
    @constraint(wm.model, q <= sum_demand)
end

function constraint_no_good{T}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    yp_solution = getvalue(wm.var[:nw][n][:yp])
    one_vars = [wm.var[:nw][n][:yp][idx[1]] for idx in keys(yp_solution) if yp_solution[idx[1]] >= 1]
    zero_vars = [wm.var[:nw][n][:yp][idx[1]] for idx in keys(yp_solution) if yp_solution[idx[1]] <= 0]
    @constraint(wm.model, sum(zero_vars) - sum(one_vars) >= 1 - length(one_vars))
end

function constraint_degree_two{T}(wm::GenericWaterModel{T}, idx, n::Int = wm.cnw)
    first = nothing
    last = nothing

    for i in wm.ref[:nw][n][:junction_connections][idx]
        connection = wm.ref[:nw][n][:connection][i]

        if connection["node1"] == idx
            other = connection["node2"]
        else
            other = connection["node1"]
        end

        if first == nothing
            first = other
        elseif first != other
            if last != nothing && last != other
                error(string("Error: adding a degree 2 constraint to a node with degree > 2: Junction ", idx))
            end

            last = other
        end
    end

    yp_first = filter(i -> wm.ref[:nw][n][:connection][i]["node1"] == first, wm.ref[:nw][n][:junction_connections][idx])
    yn_first = filter(i -> wm.ref[:nw][n][:connection][i]["node2"] == first, wm.ref[:nw][n][:junction_connections][idx])
    yp_last  = filter(i -> wm.ref[:nw][n][:connection][i]["node2"] == last,  wm.ref[:nw][n][:junction_connections][idx])
    yn_last  = filter(i -> wm.ref[:nw][n][:connection][i]["node1"] == last,  wm.ref[:nw][n][:junction_connections][idx])

    yp = wm.var[:nw][n][:yp]
    yn = wm.var[:nw][n][:yn]

    i = idx

    if !haskey(wm.con[:nw][n], :conserve_flow1)
        wm.con[:nw][n][:conserve_flow1] = Dict{String, ConstraintRef}()
        wm.con[:nw][n][:conserve_flow2] = Dict{String, ConstraintRef}()
        wm.con[:nw][n][:conserve_flow3] = Dict{String, ConstraintRef}()
        wm.con[:nw][n][:conserve_flow4] = Dict{String, ConstraintRef}()
    end

    if length(yn_first) > 0 && length(yp_last) > 0
        for i1 in yn_first
            for i2 in yp_last
                wm.con[:nw][n][:conserve_flow1][i] = @constraint(wm.model, yn[i1] == yp[i2])
                wm.con[:nw][n][:conserve_flow2][i] = @constraint(wm.model, yp[i1] == yn[i2])
                wm.con[:nw][n][:conserve_flow3][i] = @constraint(wm.model, yn[i1] + yn[i2] == 1)
                wm.con[:nw][n][:conserve_flow4][i] = @constraint(wm.model, yp[i1] + yp[i2] == 1)
            end
        end
    end

    if length(yn_first) > 0 && length(yn_last) > 0
        for i1 in yn_first
            for i2 in yn_last
                wm.con[:nw][n][:conserve_flow1][i] = @constraint(wm.model, yn[i1] == yn[i2])
                wm.con[:nw][n][:conserve_flow2][i] = @constraint(wm.model, yp[i1] == yp[i2])
                wm.con[:nw][n][:conserve_flow3][i] = @constraint(wm.model, yn[i1] + yp[i2] == 1)
                wm.con[:nw][n][:conserve_flow4][i] = @constraint(wm.model, yp[i1] + yn[i2] == 1)
            end
        end
    end

    if length(yp_first) > 0 && length(yp_last) > 0
        for i1 in yp_first
            for i2 in yp_last
                wm.con[:nw][n][:conserve_flow1][i] = @constraint(wm.model, yp[i1] == yp[i2])
                wm.con[:nw][n][:conserve_flow2][i] = @constraint(wm.model, yn[i1] == yn[i2])
                wm.con[:nw][n][:conserve_flow3][i] = @constraint(wm.model, yp[i1] + yn[i2] == 1)
                wm.con[:nw][n][:conserve_flow4][i] = @constraint(wm.model, yn[i1] + yp[i2] == 1)
            end
        end
    end

    if length(yp_first) > 0 && length(yn_last) > 0
        for i1 in yp_first
            for i2 in yn_last
                wm.con[:nw][n][:conserve_flow1][i] = @constraint(wm.model, yp[i1] == yn[i2])
                wm.con[:nw][n][:conserve_flow2][i] = @constraint(wm.model, yn[i1] == yp[i2])
                wm.con[:nw][n][:conserve_flow3][i] = @constraint(wm.model, yp[i1] + yp[i2] == 1)
                wm.con[:nw][n][:conserve_flow4][i] = @constraint(wm.model, yn[i1] + yn[i2] == 1)
            end
        end
    end
end
