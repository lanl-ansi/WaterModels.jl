#"These problem forms use binary variables to specify flow direction."
#AbstractRelaxedForm = Union{AbstractMICPForm, AbstractMILPRForm}

#"Create variables associated with the head difference for the MICP and MILP-R problems."
#function variable_absolute_head_difference(wm::GenericWaterModel{T}, n::Int = wm.cnw) where T <: AbstractRelaxedForm
#    link_ids = collect(ids(wm, :link_unknown_direction))
#    gamma_lb = Dict([(id, 0.0) for id in link_ids])
#    gamma_ub = Dict([(id, Inf) for id in link_ids])
#
#    for (id, link) in wm.ref[:nw][n][:link_unknown_direction]
#        h_i = wm.var[:nw][n][:h][link["node1"]]
#        h_j = wm.var[:nw][n][:h][link["node2"]]
#        diff_ub = getupperbound(h_i) - getlowerbound(h_j)
#        diff_lb = getlowerbound(h_i) - getupperbound(h_j)
#        gamma_ub[id] = max(abs(diff_lb), abs(diff_ub))
#    end
#
#    # Create variables that correspond to the absolute value of the head difference.
#    wm.var[:nw][n][:gamma] = @variable(wm.model, [id in link_ids],
#                                       start = gamma_ub[id],
#                                       lowerbound = gamma_lb[id],
#                                       upperbound = gamma_ub[id],
#                                       basename = "gamma_$(n)")
#end
#
##"Create variables associated with flow directions for the MICP and MILP-R problems."
##function variable_flow_direction(wm::GenericWaterModel{T}, n::Int = wm.cnw) where T <: AbstractRelaxedForm
##    # Create variables that correspond to flow moving from i to j.
##    link_ids = collect(ids(wm, :link_unknown_direction))
##    wm.var[:nw][n][:yp] = @variable(wm.model, [a in link_ids], start = 1,
##                                    category = :Bin, basename = "yp_$(n)")
##    wm.var[:nw][n][:yn] = @variable(wm.model, [a in link_ids], start = 0,
##                                    category = :Bin, basename = "yn_$(n)")
##
##    # Fix these variables if the head bounds imply they can be fixed.
##    for (id, link) in wm.ref[:nw][n][:link_unknown_direction]
##        h_i = wm.var[:nw][n][:h][parse(Int, link["node1"])]
##        h_j = wm.var[:nw][n][:h][parse(Int, link["node2"])]
##        diff_lb = getlowerbound(h_i) - getupperbound(h_j)
##        diff_ub = getupperbound(h_i) - getlowerbound(h_j)
##
##        # Use the above head difference bounds to fix upper bounds for the
##        # direction variables when applicable.
##        if (diff_lb >= 0.0)
##            setupperbound(wm.var[:nw][n][:yn][id], 0)
##            setlowerbound(wm.var[:nw][n][:yp][id], 1)
##        elseif (diff_ub <= 0.0)
##            setlowerbound(wm.var[:nw][n][:yn][id], 1)
##            setupperbound(wm.var[:nw][n][:yp][id], 0)
##        end
##    end
##end
#
#"Create variables associated with head for the MICP and MILP-R feasibility problems."
#function variable_head(wm::GenericWaterModel{T}, n::Int = wm.cnw) where T <: AbstractRelaxedForm
#    variable_head_common(wm, n) # Create h variables.
#    variable_absolute_head_difference(wm, n) # Create gamma variables.
#    variable_flow_direction(wm, n) # Create y_p and y_n variables.
#end
#
#"Create variables associated with head for the MICP and MILP-R network expansion problems."
#function variable_head_ne(wm::GenericWaterModel{T}, n::Int = wm.cnw) where T <: AbstractRelaxedForm
#    variable_head_common(wm, n) # Create h variables.
#    variable_absolute_head_difference(wm, n) # Create gamma variables.
#    variable_flow_direction(wm, n) # Create y_p and y_n variables.
#end
#
#"Create variables associated with the pipe for the MICP and MILP-R problems."
#function variable_pipe_ne(wm::GenericWaterModel{T}, n::Int = wm.cnw) where T <: AbstractRelaxedForm
#    # Create pipe variables.
#    variable_pipe_ne_common(wm, n)
#
#    # Set up required data to initialize junction variables.
#    wm.var[:nw][n][:wp] = Dict{Int, Any}()
#    wm.var[:nw][n][:wn] = Dict{Int, Any}()
#
#    # Create binary variables associated with the combination of whether a
#    # diameter is used and the direction of flow in a pipe.
#    for (id, pipe) in wm.ref[:nw][n][:ne_pipe]
#        # Obtain the diameters that can be used for the pipe.
#        diameters = [d["diameter"] for d in pipe["diameters"]]
#
#        # Instantiate the variables for the diameter/y_p combinations.
#        wm.var[:nw][n][:wp][id] = @variable(wm.model, [d in diameters],
#                                            category = :Bin,
#                                            basename = "w_p_$(n)_$(id)",
#                                            start = 0)
#
#        # Instantiate the variables for the diameter/y_n combinations.
#        wm.var[:nw][n][:wn][id] = @variable(wm.model, [d in diameters],
#                                            category = :Bin,
#                                            basename = "w_n_$(n)_$(id)",
#                                            start = 0)
#
#        # Compute the lower and upper bounds for the head difference.
#        h_i = wm.var[:nw][n][:h][pipe["node1"]]
#        h_j = wm.var[:nw][n][:h][pipe["node2"]]
#        diff_lb = getlowerbound(h_i) - getupperbound(h_j)
#        diff_ub = getupperbound(h_i) - getlowerbound(h_j)
#
#        # Use the above head difference bounds to fix upper bounds for the
#        # diameter/direction variables when applicable.
#        for diameter in diameters
#            if (diff_lb >= 0.0)
#                setupperbound(wm.var[:nw][n][:wn][id][diameter], 0)
#            elseif (diff_ub <= 0.0)
#                setupperbound(wm.var[:nw][n][:wp][id][diameter], 0)
#            end
#        end
#    end
#end
#
#function constraint_restrict_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: AbstractRelaxedForm
#    # Collect variables needed for the constraint.
#    q, h_i, h_j = get_common_variables(wm, a, n)
#
#    # Collect other variables needed for the constraint.
#    y_p = wm.var[:nw][n][:yp][a]
#    y_n = wm.var[:nw][n][:yn][a]
#
#    # Get additional data related to the variables.
#    diff_lb = getlowerbound(h_i) - getupperbound(h_j)
#    diff_ub = getupperbound(h_i) - getlowerbound(h_j)
#
#    # Get the sum of all junction demands.
#    junctions = values(wm.ref[:nw][n][:junctions])
#    sum_demand = sum(junc["demand"] for junc in junctions)
#
#    # Add the required constraints to define y_p and y_n.
#    @constraint(wm.model, (y_p - 1) * sum_demand <= q)
#    @constraint(wm.model, (1 - y_n) * sum_demand >= q)
#    @constraint(wm.model, (1 - y_p) * diff_lb <= h_i - h_j)
#    @constraint(wm.model, (1 - y_n) * diff_ub >= h_i - h_j)
#    @constraint(wm.model, y_p + y_n == 1)
#end
#
#"Constraints used to define the head difference in the MICP and MILP-R problems."
#function constraint_define_gamma(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: AbstractRelaxedForm
#    # Collect variables needed for the constraint.
#    q, h_i, h_j = get_common_variables(wm, a, n)
#
#    # Collect other variables needed for the constraint.
#    gamma = wm.var[:nw][n][:gamma][a]
#    y_p = wm.var[:nw][n][:yp][a]
#    y_n = wm.var[:nw][n][:yn][a]
#
#    # Get additional data related to the variables.
#    diff_lb = getlowerbound(h_i) - getupperbound(h_j)
#    diff_ub = getupperbound(h_i) - getlowerbound(h_j)
#
#    # Add the required constraints to define gamma.
#    @constraint(wm.model, h_j - h_i + diff_lb * (y_p - y_n + 1) <= gamma)
#    @constraint(wm.model, h_i - h_j + diff_ub * (y_p - y_n - 1) <= gamma)
#    @constraint(wm.model, h_j - h_i + diff_ub * (y_p - y_n + 1) >= gamma)
#    @constraint(wm.model, h_i - h_j + diff_lb * (y_p - y_n - 1) >= gamma)
#end
#
#"Constraints used to define the head difference in the MICP and MILP-R expansion planning problems."
#function constraint_define_gamma_hw_ne(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: AbstractRelaxedForm
#    # Collect variables and parameters needed for the constraint.
#    q, h_i, h_j = get_common_variables(wm, a, n)
#    y_p = wm.var[:nw][n][:yp][a]
#    y_n = wm.var[:nw][n][:yn][a]
#
#    # Get the various pipe diameters from which we can select.
#    pipe = wm.ref[:nw][n][:ne_pipe][a]
#    diameters = [d["diameter"] for d in pipe["diameters"]]
#
#    # Add a constraint that says at most one diameter must be selected.
#    @constraint(wm.model, sum(wm.var[:nw][n][:psi][a][d] for d in diameters) == 1)
#
#    # Get the sum of all junction demands.
#    junctions = values(wm.ref[:nw][n][:junctions])
#    sum_demand = sum(junction["demand"] for junction in junctions)
#    @constraint(wm.model, sum(wm.var[:nw][n][:wp][a][d] for d in diameters) == y_p)
#    @constraint(wm.model, sum(wm.var[:nw][n][:wn][a][d] for d in diameters) == y_n)
#    @constraint(wm.model, y_p + y_n == 1)
#    @constraint(wm.model, (y_p - 1) * sum_demand <= q)
#    @constraint(wm.model, (1 - y_n) * sum_demand >= q)
#
#    # Get the pipe associated with the pipe index a.
#    pipe = wm.ref[:nw][n][:ne_pipe][a]
#
#    # Constrain each gamma variable.
#    for diameter in diameters
#        # Gather the required variables and constants.
#        psi_d = wm.var[:nw][n][:psi][a][diameter]
#        gamma_d = wm.var[:nw][n][:gamma][a][diameter]
#        z_p = wm.var[:nw][n][:wp][a][diameter]
#        z_n = wm.var[:nw][n][:wn][a][diameter]
#        lambda_d = calc_friction_factor_hw_ne(pipe, diameter)
#        gamma_ub = (getupperbound(h_i) - getlowerbound(h_j)) / lambda_d
#        gamma_lb = (getlowerbound(h_i) - getupperbound(h_j)) / lambda_d
#
#        # Ensure that z is chosen when the diameter is chosen.
#        @constraint(wm.model, (z_p + z_n) == psi_d)
#
#        # Link z_p and z_n with q and (h_i - h_j).
#        @constraint(wm.model, (z_p - 1) * sum_demand <= q)
#        @constraint(wm.model, (1 - z_n) * sum_demand >= q)
#        @constraint(wm.model, (1 - z_p) * gamma_lb <= (h_i - h_j) / lambda_d)
#        @constraint(wm.model, (1 - z_n) * gamma_ub >= (h_i - h_j) / lambda_d)
#        @constraint(wm.model, (z_p + z_n) * gamma_ub >= gamma_d)
#        @constraint(wm.model, (z_p + z_n) * gamma_lb <= gamma_d)
#
#        # Create McCormick constraints that ensure gamma_d is positive.
#        @constraint(wm.model, (h_j - h_i) / lambda_d + gamma_lb * (z_p - z_n + 1) <= gamma_d)
#        @constraint(wm.model, (h_i - h_j) / lambda_d + gamma_ub * (z_p - z_n - 1) <= gamma_d)
#        @constraint(wm.model, (h_j - h_i) / lambda_d + gamma_ub * (z_p - z_n + 1) >= gamma_d)
#        @constraint(wm.model, (h_i - h_j) / lambda_d + gamma_lb * (z_p - z_n - 1) >= gamma_d)
#    end
#end
#
#function set_initial_solution_ne(wm::GenericWaterModel{T}, wm_solved::GenericWaterModel) where T <: AbstractRelaxedForm
#    for i in [collect(ids(wm, :junctions)); collect(ids(wm, :reservoirs))]
#        h_i_sol = getvalue(wm_solved.var[:nw][wm_solved.cnw][:h][i])
#        setvalue(wm.var[:nw][wm.cnw][:h][i], h_i_sol)
#    end
#
#    objective_value = 0.0
#
#    for ij in collect(ids(wm, :link))
#        i = parse(Int, wm_solved.ref[:nw][wm_solved.cnw][:pipes][ij]["node1"])
#        j = parse(Int, wm_solved.ref[:nw][wm_solved.cnw][:pipes][ij]["node2"])
#        q_ij_sol = getvalue(wm_solved.var[:nw][wm_solved.cnw][:q][ij])
#        h_i_sol = getvalue(wm_solved.var[:nw][wm_solved.cnw][:h][i])
#        h_j_sol = getvalue(wm_solved.var[:nw][wm_solved.cnw][:h][j])
#        setvalue(wm.var[:nw][wm.cnw][:q][ij], q_ij_sol)
#
#        flow_direction = sign(q_ij_sol)
#        setvalue(wm.var[:nw][wm.cnw][:yp][ij], 1 * (flow_direction == 1))
#        setvalue(wm.var[:nw][wm.cnw][:yn][ij], 1 * (flow_direction == -1))
#
#        # Get the various pipe diameters from which we can select.
#        pipe = wm.ref[:nw][wm.cnw][:ne_pipe][ij]
#        diameters = [d["diameter"] for d in pipe["diameters"]]
#        diameter_sol = wm_solved.ref[:nw][wm_solved.cnw][:pipes][ij]["diameter"]
#
#        for diameter in diameters
#            if diameter == diameter_sol
#                pipe = wm.ref[:nw][wm.cnw][:ne_pipe][ij]
#                pipe_sol = wm_solved.ref[:nw][wm_solved.cnw][:pipes][ij]
#                lambda = calc_friction_factor_hw_ne(pipe_sol, diameter)
#                gamma = flow_direction * (h_i_sol - h_j_sol) / lambda
#                cost_per_unit_length = [d["costPerUnitLength"] for d in filter((d) -> d["diameter"] == diameter_sol, pipe["diameters"])][1]
#                objective_value += pipe["length"] * cost_per_unit_length
#                setvalue(wm.var[:nw][wm.cnw][:gamma][ij][diameter], gamma)
#                setvalue(wm.var[:nw][wm.cnw][:gamma_sum][ij], gamma)
#                setvalue(wm.var[:nw][wm.cnw][:psi][ij][diameter], 1)
#                setvalue(wm.var[:nw][wm.cnw][:wp][ij][diameter], 1 * (flow_direction == 1))
#                setvalue(wm.var[:nw][wm.cnw][:wn][ij][diameter], 1 * (flow_direction == -1))
#            else
#                setvalue(wm.var[:nw][wm.cnw][:gamma][ij][diameter], 0.0)
#                setvalue(wm.var[:nw][wm.cnw][:psi][ij][diameter], 0)
#                setvalue(wm.var[:nw][wm.cnw][:wp][ij][diameter], 0)
#                setvalue(wm.var[:nw][wm.cnw][:wn][ij][diameter], 0)
#            end
#        end
#    end
#
#    setvalue(wm.var[:nw][wm.cnw][:objective], objective_value * 1.0e-6)
#end
#
#function set_initial_solution_from_cvx(wm::GenericWaterModel{T}, cvx::GenericWaterModel) where T <: AbstractRelaxedForm
#    h_sol = Dict{Int, Float64}()
#
#    num_junctions = length(cvx.ref[:nw][cvx.cnw][:junctions])
#    num_reservoirs = length(cvx.ref[:nw][cvx.cnw][:reservoirs])
#    num_nodes = num_junctions + num_reservoirs
#    num_arcs = length(cvx.ref[:nw][cvx.cnw][:links])
#    A = zeros(Float64, num_arcs + num_reservoirs, num_nodes)
#    b = zeros(Float64, num_arcs + num_reservoirs, 1)
#
#    for (a, link) in cvx.ref[:nw][cvx.cnw][:links]
#        A[a, link["node1"]] = 1
#        A[a, link["node2"]] = -1
#        q_p = getvalue(cvx.var[:nw][cvx.cnw][:qp][a])
#        q_n = getvalue(cvx.var[:nw][cvx.cnw][:qn][a])
#        r = calc_resistance_per_length_hw(link)
#        b[a] = link["length"] * r * (q_p - q_n)
#    end
#
#    k = 1
#    for (i, reservoir) in cvx.ref[:nw][cvx.cnw][:reservoirs]
#        row = num_arcs + k
#        A[row, i] = 1
#        b[row] = reservoir["head"]
#    end
#
#    h_sol = A \ b
#
#    for (i, junction) in wm.ref[:nw][wm.cnw][:junctions]
#        #h_sol[i] = getdual(cvx.con[:nw][cvx.cnw][:flow_conservation][i])
#        set_start_value(wm.var[:nw][wm.cnw][:h][i], h_sol[i])
#    end
#
#    for (i, reservoir) in wm.ref[:nw][wm.cnw][:reservoirs]
#        #h_sol[i] = reservoir["head"]
#        set_start_value(wm.var[:nw][wm.cnw][:h][i], h_sol[i])
#    end
#
#    q_sol = Dict{Int, Float64}()
#
#    for (ij, link) in wm.ref[:nw][wm.cnw][:links]
#        i = parse(Int, cvx.ref[:nw][cvx.cnw][:pipes][ij]["node1"])
#        j = parse(Int, cvx.ref[:nw][cvx.cnw][:pipes][ij]["node2"])
#        q_p_ij = getvalue(cvx.var[:nw][cvx.cnw][:qp][ij])
#        q_n_ij = getvalue(cvx.var[:nw][cvx.cnw][:qn][ij])
#        q_sol[ij] = q_p_ij - q_n_ij
#        L = link["length"]
#        r = calc_resistance_per_length_hw(link)
#    end
#
#    objective_value = 0.0
#
#    for ij in collect(ids(wm, :link))
#        i = parse(Int, cvx.ref[:nw][cvx.cnw][:pipes][ij]["node1"])
#        j = parse(Int, cvx.ref[:nw][cvx.cnw][:pipes][ij]["node2"])
#
#        flow_direction = sign(q_sol[ij])
#        set_start_value(wm.var[:nw][wm.cnw][:q][ij], q_sol[ij])
#        set_start_value(wm.var[:nw][wm.cnw][:yp][ij], 1 * (flow_direction == 1))
#        set_start_value(wm.var[:nw][wm.cnw][:yn][ij], 1 * (flow_direction == -1))
#
#        # Get the various pipe diameters from which we can select.
#        pipe = wm.ref[:nw][wm.cnw][:ne_pipe][ij]
#        diameters = [d["diameter"] for d in pipe["diameters"]]
#        diameter_sol = cvx.ref[:nw][cvx.cnw][:pipes][ij]["diameter"]
#
#        for diameter in diameters
#            if diameter == diameter_sol
#                pipe = wm.ref[:nw][wm.cnw][:ne_pipe][ij]
#                pipe_sol = cvx.ref[:nw][cvx.cnw][:pipes][ij]
#                lambda = calc_friction_factor_hw_ne(pipe_sol, diameter)
#                gamma = flow_direction * (h_sol[i] - h_sol[j]) / lambda
#                cost_per_unit_length = [d["costPerUnitLength"] for d in filter((d) -> d["diameter"] == diameter_sol, pipe["diameters"])][1]
#                objective_value += pipe["length"] * cost_per_unit_length
#                set_start_value(wm.var[:nw][wm.cnw][:gamma][ij][diameter], gamma)
#                set_start_value(wm.var[:nw][wm.cnw][:gamma_sum][ij], gamma)
#                set_start_value(wm.var[:nw][wm.cnw][:psi][ij][diameter], 1)
#                set_start_value(wm.var[:nw][wm.cnw][:wp][ij][diameter], 1 * (flow_direction == 1))
#                set_start_value(wm.var[:nw][wm.cnw][:wn][ij][diameter], 1 * (flow_direction == -1))
#            else
#                set_start_value(wm.var[:nw][wm.cnw][:gamma][ij][diameter], 0.0)
#                set_start_value(wm.var[:nw][wm.cnw][:psi][ij][diameter], 0)
#                set_start_value(wm.var[:nw][wm.cnw][:wp][ij][diameter], 0)
#                set_start_value(wm.var[:nw][wm.cnw][:wn][ij][diameter], 0)
#            end
#        end
#    end
#
#    set_start_value(wm.var[:nw][wm.cnw][:objective], objective_value * 1.0e-6)
#end
#
#"Constraint to ensure at least one direction is set to take flow away from a source."
#function constraint_source_flow(wm::GenericWaterModel{T}, i::Int, f_branches, t_branches, n::Int = wm.cnw) where T <: AbstractRelaxedForm
#    dirs_out_1 = Array{JuMP.Variable}([wm.var[:nw][n][:yp][a] for a in f_branches])
#    dirs_out_2 = Array{JuMP.Variable}([wm.var[:nw][n][:yn][a] for a in t_branches])
#    @constraint(wm.model, sum(dirs_out_1) + sum(dirs_out_2) >= 1)
#end
#
#"Constraint to ensure at least one direction is set to take flow to a junction with demand."
#function constraint_sink_flow(wm::GenericWaterModel{T}, i::Int, f_branches, t_branches, n::Int = wm.cnw) where T <: AbstractRelaxedForm
#    dirs_in_1 = Array{JuMP.Variable}([wm.var[:nw][n][:yp][a] for a in t_branches])
#    dirs_in_2 = Array{JuMP.Variable}([wm.var[:nw][n][:yn][a] for a in f_branches])
#    @constraint(wm.model, sum(dirs_in_1) + sum(dirs_in_2) >= 1)
#end
#
#function constraint_junction_mass_flow(wm::GenericWaterModel{T}, i::Int, n::Int = wm.cnw) where T <: AbstractRelaxedForm
#    constraint_flow_conservation(wm, i, n)
#
#    is_reservoir = haskey(wm.ref[:nw][n][:reservoirs], i)
#    is_junction = haskey(wm.ref[:nw][n][:junctions], i)
#
#    has_demand = false
#    if is_junction
#        has_demand = wm.ref[:nw][n][:junctions][i]["demand"] > 0.0
#    end
#
#    if is_reservoir && !has_demand
#        constraint_source_flow(wm, i, n)
#    elseif !is_reservoir && has_demand
#        constraint_sink_flow(wm, i, n)
#    end
#
#    node_degree = length(wm.ref[:nw][n][:junction_links][i])
#    if !is_reservoir && !has_demand && node_degree == 2
#        constraint_degree_two(wm, i, n)
#    end
#end
#
#function constraint_no_good_ne(wm::GenericWaterModel{T}, n::Int = wm.cnw) where T <: AbstractRelaxedForm
#    pipe_ids = ids(wm, :ne_pipe)
#    psi_vars = Array{JuMP.Variable}([wm.var[:nw][n][:psi][a][d["diameter"]] for a in ids(wm, :ne_pipe) for d in wm.ref[:nw][n][:ne_pipe][a]["diameters"]])
#    psi_ones = Array{JuMP.Variable}([psi for psi in psi_vars if getvalue(psi) > 1.0e-4])
#    psi_zeros = Array{JuMP.Variable}([psi for psi in psi_vars if getvalue(psi) <= 1.0e-4])
#    @constraint(wm.model, sum(psi_zeros) - sum(psi_ones) >= 1 - length(psi_ones))
#end
#
#function constraint_no_good_ne_cb(wm::GenericWaterModel{T}, cb, n::Int = wm.cnw) where T <: AbstractRelaxedForm
#    pipe_ids = ids(wm, :ne_pipe)
#    psi_vars = Array{JuMP.Variable}([wm.var[:nw][n][:psi][a][d["diameter"]] for a in ids(wm, :ne_pipe) for d in wm.ref[:nw][n][:ne_pipe][a]["diameters"]])
#    psi_ones = Array{JuMP.Variable}([psi for psi in psi_vars if getvalue(psi) > 1.0e-4])
#    psi_zeros = Array{JuMP.Variable}([psi for psi in psi_vars if getvalue(psi) <= 1.0e-4])
#    @lazyconstraint(cb, sum(psi_zeros) - sum(psi_ones) >= 1 - length(psi_ones))
#end
#
#function constraint_degree_two(wm::GenericWaterModel{T}, idx, n::Int = wm.cnw) where T <: AbstractRelaxedForm
#    first = nothing
#    last = nothing
#
#    for i in wm.ref[:nw][n][:junction_links][idx]
#        link = wm.ref[:nw][n][:links][i]
#
#        if parse(Int, link["node1"]) == idx
#            other = link["node2"]
#        else
#            other = link["node1"]
#        end
#
#        if first == nothing
#            first = other
#        elseif first != other
#            if last != nothing && last != other
#                error(string("Error: adding a degree 2 constraint to a node with degree > 2: Junction ", idx))
#            end
#
#            last = other
#        end
#    end
#
#    yp_first = filter(i -> wm.ref[:nw][n][:links][i]["node1"] == first, wm.ref[:nw][n][:junction_links][idx])
#    yn_first = filter(i -> wm.ref[:nw][n][:links][i]["node2"] == first, wm.ref[:nw][n][:junction_links][idx])
#    yp_last  = filter(i -> wm.ref[:nw][n][:links][i]["node2"] == last,  wm.ref[:nw][n][:junction_links][idx])
#    yn_last  = filter(i -> wm.ref[:nw][n][:links][i]["node1"] == last,  wm.ref[:nw][n][:junction_links][idx])
#
#    yp = wm.var[:nw][n][:yp]
#    yn = wm.var[:nw][n][:yn]
#
#    i = idx
#
#    if !haskey(wm.con[:nw][n], :conserve_flow1)
#        wm.con[:nw][n][:conserve_flow1] = Dict{Int, ConstraintRef}()
#        wm.con[:nw][n][:conserve_flow2] = Dict{Int, ConstraintRef}()
#        wm.con[:nw][n][:conserve_flow3] = Dict{Int, ConstraintRef}()
#        wm.con[:nw][n][:conserve_flow4] = Dict{Int, ConstraintRef}()
#    end
#
#    if length(yn_first) > 0 && length(yp_last) > 0
#        for i1 in yn_first
#            for i2 in yp_last
#                wm.con[:nw][n][:conserve_flow1][i] = @constraint(wm.model, yn[i1] == yp[i2])
#                wm.con[:nw][n][:conserve_flow2][i] = @constraint(wm.model, yp[i1] == yn[i2])
#                wm.con[:nw][n][:conserve_flow3][i] = @constraint(wm.model, yn[i1] + yn[i2] == 1)
#                wm.con[:nw][n][:conserve_flow4][i] = @constraint(wm.model, yp[i1] + yp[i2] == 1)
#            end
#        end
#    end
#
#    if length(yn_first) > 0 && length(yn_last) > 0
#        for i1 in yn_first
#            for i2 in yn_last
#                wm.con[:nw][n][:conserve_flow1][i] = @constraint(wm.model, yn[i1] == yn[i2])
#                wm.con[:nw][n][:conserve_flow2][i] = @constraint(wm.model, yp[i1] == yp[i2])
#                wm.con[:nw][n][:conserve_flow3][i] = @constraint(wm.model, yn[i1] + yp[i2] == 1)
#                wm.con[:nw][n][:conserve_flow4][i] = @constraint(wm.model, yp[i1] + yn[i2] == 1)
#            end
#        end
#    end
#
#    if length(yp_first) > 0 && length(yp_last) > 0
#        for i1 in yp_first
#            for i2 in yp_last
#                wm.con[:nw][n][:conserve_flow1][i] = @constraint(wm.model, yp[i1] == yp[i2])
#                wm.con[:nw][n][:conserve_flow2][i] = @constraint(wm.model, yn[i1] == yn[i2])
#                wm.con[:nw][n][:conserve_flow3][i] = @constraint(wm.model, yp[i1] + yn[i2] == 1)
#                wm.con[:nw][n][:conserve_flow4][i] = @constraint(wm.model, yn[i1] + yp[i2] == 1)
#            end
#        end
#    end
#
#    if length(yp_first) > 0 && length(yn_last) > 0
#        for i1 in yp_first
#            for i2 in yn_last
#                wm.con[:nw][n][:conserve_flow1][i] = @constraint(wm.model, yp[i1] == yn[i2])
#                wm.con[:nw][n][:conserve_flow2][i] = @constraint(wm.model, yn[i1] == yp[i2])
#                wm.con[:nw][n][:conserve_flow3][i] = @constraint(wm.model, yp[i1] + yp[i2] == 1)
#                wm.con[:nw][n][:conserve_flow4][i] = @constraint(wm.model, yn[i1] + yn[i2] == 1)
#            end
#        end
#    end
#end
