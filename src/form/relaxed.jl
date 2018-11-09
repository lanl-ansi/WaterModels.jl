"These problem forms use binary variables to specify flow direction."
AbstractRelaxedForm = Union{AbstractMICPForm, AbstractMILPRForm}

"Create variables associated with the head difference for the MICP and MILP-R problems."
function variable_absolute_head_difference{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    connections = wm.ref[:nw][n][:connection_unknown_direction]
    gamma_lb = Dict([(id, 0.0) for id in keys(connections)])
    gamma_ub = Dict([(id, Inf) for id in keys(connections)])

    for (id, connection) in connections
        h_i = wm.var[:nw][n][:h][connection["node1"]]
        h_j = wm.var[:nw][n][:h][connection["node2"]]
        diff_ub = getupperbound(h_i) - getlowerbound(h_j)
        diff_lb = getlowerbound(h_i) - getupperbound(h_j)
        gamma_ub[id] = max(abs(diff_lb), abs(diff_ub))
    end

    # Create variables that correspond to the absolute value of the head difference.
    wm.var[:nw][n][:gamma] = @variable(wm.model, [id in keys(connections)],
                                       lowerbound = gamma_lb[id],
                                       upperbound = gamma_ub[id],
                                       basename = "gamma_$(n)")
end

"Create variables associated with flow directions for the MICP and MILP-R problems."
function variable_flow_direction{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    connections = wm.ref[:nw][n][:connection_unknown_direction]

    # Create variables that correspond to flow moving from i to j.
    wm.var[:nw][n][:yp] = @variable(wm.model, [id in keys(connections)],
                                    category = :Bin, basename = "yp_$(n)")
    wm.var[:nw][n][:yn] = @variable(wm.model, [id in keys(connections)],
                                    category = :Bin, basename = "yp_$(n)")

    # Fix these variables if the head bounds imply they can be fixed.
    for (id, connection) in connections
        h_i = wm.var[:nw][n][:h][connection["node1"]]
        h_j = wm.var[:nw][n][:h][connection["node2"]]
        diff_lb = getlowerbound(h_i) - getupperbound(h_j)
        diff_ub = getupperbound(h_i) - getlowerbound(h_j)

        # Use the above head difference bounds to fix upper bounds for the
        # direction variables when applicable.
        if (diff_lb >= 0.0)
            setupperbound(wm.var[:nw][n][:yn][id], 0)
            setlowerbound(wm.var[:nw][n][:yp][id], 1)
        elseif (diff_ub <= 0.0)
            setlowerbound(wm.var[:nw][n][:yn][id], 1)
            setupperbound(wm.var[:nw][n][:yp][id], 0)
        end
    end
end

"Create variables associated with head for the MICP and MILP-R feasibility problems."
function variable_head{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    variable_head_common(wm, n) # Create h variables.
    variable_absolute_head_difference(wm, n) # Create gamma variables.
    variable_flow_direction(wm, n) # Create y_p and y_n variables.
end

"Create variables associated with head for the MICP and MILP-R network expansion problems."
function variable_head_ne{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    variable_head_common(wm, n) # Create h variables.
    variable_absolute_head_difference(wm, n) # Create gamma variables.
    variable_flow_direction(wm, n) # Create y_p and y_n variables.
end

"Create variables associated with the pipe for the MICP and MILP-R problems."
function variable_pipe_ne{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    # Create pipe variables.
    variable_pipe_ne_common(wm, n)

    # Set up required data to initialize junction variables.
    wm.var[:nw][n][:wp] = Dict{String, Any}()
    wm.var[:nw][n][:wn] = Dict{String, Any}()

    # Create binary variables associated with the combination of whether a
    # diameter is used and the direction of flow in a pipe.
    for (id, pipe) in wm.ref[:nw][n][:ne_pipe]
        # Obtain the diameters that can be used for the pipe.
        diameters = [d["diameter"] for d in pipe["diameters"]]

        # Instantiate the variables for the diameter/y_p combinations.
        wm.var[:nw][n][:wp][id] = @variable(wm.model, [d in diameters],
                                            category = :Bin,
                                            basename = "w_p_$(n)_$(id)",
                                            start = 0)

        # Instantiate the variables for the diameter/y_n combinations.
        wm.var[:nw][n][:wn][id] = @variable(wm.model, [d in diameters],
                                            category = :Bin,
                                            basename = "w_n_$(n)_$(id)",
                                            start = 0)

        # Compute the lower and upper bounds for the head difference.
        h_i = wm.var[:nw][n][:h][pipe["node1"]]
        h_j = wm.var[:nw][n][:h][pipe["node2"]]
        diff_lb = getlowerbound(h_i) - getupperbound(h_j)
        diff_ub = getupperbound(h_i) - getlowerbound(h_j)

        # Use the above head difference bounds to fix upper bounds for the
        # diameter/direction variables when applicable.
        for diameter in diameters
            if (diff_lb >= 0.0)
                setupperbound(wm.var[:nw][n][:wn][id][diameter], 0)
            elseif (diff_ub <= 0.0)
                setupperbound(wm.var[:nw][n][:wp][id][diameter], 0)
            end
        end
    end
end

function constraint_restrict_direction{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

    # Collect other variables needed for the constraint.
    y_p = wm.var[:nw][n][:yp][a]
    y_n = wm.var[:nw][n][:yn][a]

    # Get additional data related to the variables.
    diff_lb = getlowerbound(h_i) - getupperbound(h_j)
    diff_ub = getupperbound(h_i) - getlowerbound(h_j)

    # Get the sum of all junction demands.
    junctions = values(wm.ref[:nw][n][:junctions])
    sum_demand = sum(junc["demand"] for junc in junctions)

    # Add the required constraints to define y_p and y_n.
    @constraint(wm.model, (y_p - 1) * sum_demand <= q)
    @constraint(wm.model, (1 - y_n) * sum_demand >= q)
    @constraint(wm.model, (1 - y_p) * diff_lb <= h_i - h_j)
    @constraint(wm.model, (1 - y_n) * diff_ub >= h_i - h_j)
    @constraint(wm.model, y_p + y_n == 1)
end

"Constraints used to define the head difference in the MICP and MILP-R problems."
function constraint_define_gamma{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

    # Collect other variables needed for the constraint.
    gamma = wm.var[:nw][n][:gamma][a]
    y_p = wm.var[:nw][n][:yp][a]
    y_n = wm.var[:nw][n][:yn][a]

    # Get additional data related to the variables.
    diff_lb = getlowerbound(h_i) - getupperbound(h_j)
    diff_ub = getupperbound(h_i) - getlowerbound(h_j)

    # Add the required constraints to define gamma.
    @constraint(wm.model, h_j - h_i + diff_lb * (y_p - y_n + 1) <= gamma)
    @constraint(wm.model, h_i - h_j + diff_ub * (y_p - y_n - 1) <= gamma)
    @constraint(wm.model, h_j - h_i + diff_ub * (y_p - y_n + 1) >= gamma)
    @constraint(wm.model, h_i - h_j + diff_lb * (y_p - y_n - 1) >= gamma)
end

"Constraints used to define the head difference in the MICP and MILP-R expansion planning problems."
function constraint_define_gamma_hw_ne{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)
    y_p = wm.var[:nw][n][:yp][a]
    y_n = wm.var[:nw][n][:yn][a]

    # Add a constraint that says at most one diameter must be selected.
    @constraint(wm.model, sum(wm.var[:nw][n][:psi][a]) == 1)

    # Get the sum of all junction demands.
    junctions = values(wm.ref[:nw][n][:junctions])
    sum_demand = sum(junction["demand"] for junction in junctions)
    @constraint(wm.model, sum(wm.var[:nw][n][:wp][a]) == y_p)
    @constraint(wm.model, sum(wm.var[:nw][n][:wn][a]) == y_n)
    @constraint(wm.model, y_p + y_n == 1)
    @constraint(wm.model, (y_p - 1) * sum_demand <= q)
    @constraint(wm.model, (1 - y_n) * sum_demand >= q)

    # Get the pipe associated with the pipe index a.
    pipe = wm.ref[:nw][n][:ne_pipe][a]

    # Get the various pipe diameters from which we can select.
    diameters = [key[1] for key in keys(wm.var[:nw][n][:psi][a])]

    # Constrain each gamma variable.
    for diameter in diameters
        # Gather the required variables and constants.
        psi_d = wm.var[:nw][n][:psi][a][diameter]
        gamma_d = wm.var[:nw][n][:gamma][a][diameter]
        z_p = wm.var[:nw][n][:wp][a][diameter]
        z_n = wm.var[:nw][n][:wn][a][diameter]
        lambda_d = calc_friction_factor_hw_ne(pipe, diameter)
        gamma_ub = (getupperbound(h_i) - getlowerbound(h_j)) / lambda_d
        gamma_lb = (getlowerbound(h_i) - getupperbound(h_j)) / lambda_d

        # Ensure that z is chosen when the diameter is chosen.
        @constraint(wm.model, (z_p + z_n) == psi_d)

        # Link z_p and z_n with q and (h_i - h_j).
        @constraint(wm.model, (z_p - 1) * sum_demand <= q)
        @constraint(wm.model, (1 - z_n) * sum_demand >= q)
        @constraint(wm.model, (1 - z_p) * gamma_lb <= (h_i - h_j) / lambda_d)
        @constraint(wm.model, (1 - z_n) * gamma_ub >= (h_i - h_j) / lambda_d)
        @constraint(wm.model, (z_p + z_n) * gamma_ub >= gamma_d)
        @constraint(wm.model, (z_p + z_n) * gamma_lb <= gamma_d)

        # Create McCormick constraints that ensure gamma_d is positive.
        @constraint(wm.model, (h_j - h_i) / lambda_d + gamma_lb * (z_p - z_n + 1) <= gamma_d)
        @constraint(wm.model, (h_i - h_j) / lambda_d + gamma_ub * (z_p - z_n - 1) <= gamma_d)
        @constraint(wm.model, (h_j - h_i) / lambda_d + gamma_ub * (z_p - z_n + 1) >= gamma_d)
        @constraint(wm.model, (h_i - h_j) / lambda_d + gamma_lb * (z_p - z_n - 1) >= gamma_d)
    end
end

function set_initial_solution_ne{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, wm_solved::GenericWaterModel)
    for i in [key[1] for key in keys(wm_solved.var[:nw][wm_solved.cnw][:h])]
        h_i_sol = getvalue(wm_solved.var[:nw][wm_solved.cnw][:h][i])
        setvalue(wm.var[:nw][wm.cnw][:h][i], h_i_sol)
    end

    objective_value = 0.0

    for ij in [key[1] for key in keys(wm_solved.var[:nw][wm_solved.cnw][:q])]
        i = wm_solved.ref[:nw][wm_solved.cnw][:pipes][ij]["node1"]
        j = wm_solved.ref[:nw][wm_solved.cnw][:pipes][ij]["node2"]
        q_ij_sol = getvalue(wm_solved.var[:nw][wm_solved.cnw][:q][ij])
        h_i_sol = getvalue(wm_solved.var[:nw][wm_solved.cnw][:h][i])
        h_j_sol = getvalue(wm_solved.var[:nw][wm_solved.cnw][:h][j])
        setvalue(wm.var[:nw][wm.cnw][:q][ij], q_ij_sol)

        flow_direction = sign(q_ij_sol)
        setvalue(wm.var[:nw][wm.cnw][:yp][ij], 1 * (flow_direction == 1))
        setvalue(wm.var[:nw][wm.cnw][:yn][ij], 1 * (flow_direction == -1))

        diameter_sol = wm_solved.ref[:nw][wm_solved.cnw][:pipes][ij]["diameter"]
        diameters = [key[1] for key in keys(wm.var[:nw][wm.cnw][:psi][ij])]

        for diameter in diameters
            if diameter == diameter_sol
                pipe = wm.ref[:nw][wm.cnw][:ne_pipe][ij]
                pipe_sol = wm_solved.ref[:nw][wm_solved.cnw][:pipes][ij]
                lambda = calc_friction_factor_hw_ne(pipe_sol, diameter)
                gamma = flow_direction * (h_i_sol - h_j_sol) / lambda
                cost_per_unit_length = [d["costPerUnitLength"] for d in filter((d) -> d["diameter"] == diameter_sol, pipe["diameters"])][1]
                objective_value += pipe["length"] * cost_per_unit_length
                setvalue(wm.var[:nw][wm.cnw][:gamma][ij][diameter], gamma)
                setvalue(wm.var[:nw][wm.cnw][:gamma_sum][ij], gamma)
                setvalue(wm.var[:nw][wm.cnw][:psi][ij][diameter], 1)
                setvalue(wm.var[:nw][wm.cnw][:wp][ij][diameter], 1 * (flow_direction == 1))
                setvalue(wm.var[:nw][wm.cnw][:wn][ij][diameter], 1 * (flow_direction == -1))
            else
                setvalue(wm.var[:nw][wm.cnw][:gamma][ij][diameter], 0.0)
                setvalue(wm.var[:nw][wm.cnw][:psi][ij][diameter], 0)
                setvalue(wm.var[:nw][wm.cnw][:wp][ij][diameter], 0)
                setvalue(wm.var[:nw][wm.cnw][:wn][ij][diameter], 0)
            end
        end
    end

    setvalue(wm.var[:nw][wm.cnw][:objective], objective_value * 1.0e-6)
end

"Constraint to ensure at least one direction is set to take flow away from a source."
function constraint_source_flow{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, i, f_branches, t_branches, n::Int = wm.cnw)
    dirs_out_1 = Array{JuMP.Variable}([wm.var[:nw][n][:yp][a] for a in f_branches])
    dirs_out_2 = Array{JuMP.Variable}([wm.var[:nw][n][:yn][a] for a in t_branches])
    @constraint(wm.model, sum(dirs_out_1) + sum(dirs_out_2) >= 1)
end

"Constraint to ensure at least one direction is set to take flow to a junction with demand."
function constraint_sink_flow{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, i, f_branches, t_branches, n::Int = wm.cnw)
    dirs_in_1 = Array{JuMP.Variable}([wm.var[:nw][n][:yp][a] for a in t_branches])
    dirs_in_2 = Array{JuMP.Variable}([wm.var[:nw][n][:yn][a] for a in f_branches])
    @constraint(wm.model, sum(dirs_in_1) + sum(dirs_in_2) >= 1)
end

function constraint_junction_mass_flow{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, i, n::Int = wm.cnw)
    constraint_flow_conservation(wm, i, n)

    is_reservoir = haskey(wm.ref[:nw][n][:reservoirs], i)
    is_junction = haskey(wm.ref[:nw][n][:junctions], i)

    has_demand = false
    if is_junction
        has_demand = wm.ref[:nw][n][:junctions][i]["demand"] > 0.0
    end

    if is_reservoir && !has_demand
        constraint_source_flow(wm, i, n)
    elseif !is_reservoir && has_demand
        constraint_sink_flow(wm, i, n)
    end

    node_degree = length(wm.ref[:nw][n][:junction_connections][i])
    if !is_reservoir && !has_demand && node_degree == 2
        constraint_degree_two(wm, i, n)
    end
end

function constraint_no_good_ne{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    pipe_ids = ids(wm, :ne_pipe)
    psi_vars = Array{JuMP.Variable}([wm.var[:nw][n][:psi][a][d["diameter"]] for a in ids(wm, :ne_pipe) for d in wm.ref[:nw][n][:ne_pipe][a]["diameters"]])
    psi_ones = Array{JuMP.Variable}([psi for psi in psi_vars if getvalue(psi) > 1.0e-4])
    psi_zeros = Array{JuMP.Variable}([psi for psi in psi_vars if getvalue(psi) <= 1.0e-4])
    @constraint(wm.model, sum(psi_zeros) - sum(psi_ones) >= 1 - length(psi_ones))
end

function constraint_no_good_ne_cb{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, cb, n::Int = wm.cnw)
    pipe_ids = ids(wm, :ne_pipe)
    psi_vars = Array{JuMP.Variable}([wm.var[:nw][n][:psi][a][d["diameter"]] for a in ids(wm, :ne_pipe) for d in wm.ref[:nw][n][:ne_pipe][a]["diameters"]])
    psi_ones = Array{JuMP.Variable}([psi for psi in psi_vars if getvalue(psi) > 1.0e-4])
    psi_zeros = Array{JuMP.Variable}([psi for psi in psi_vars if getvalue(psi) <= 1.0e-4])
    @lazyconstraint(cb, sum(psi_zeros) - sum(psi_ones) >= 1 - length(psi_ones))
end

function constraint_degree_two{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, idx, n::Int = wm.cnw)
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
