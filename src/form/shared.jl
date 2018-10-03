"Set new bounds for q given some specified direction of flow (-1 or 1)."
function fix_flow_direction(q::JuMP.Variable, direction::Int)
    # Fix the direction of the flow.
    setlowerbound(q, direction == 1 ? 0.0 : getlowerbound(q))
    setupperbound(q, direction == 1 ? getupperbound(q) : 0.0)
end

"Get variables commonly used in the construction of head loss constraints."
function get_common_variables{T <: AbstractWaterFormulation}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Get source and target nodes corresponding to the edge.
    i = wm.ref[:nw][n][:pipes][a]["node1"]
    j = wm.ref[:nw][n][:pipes][a]["node2"]

    # Collect variables needed for the constraint.
    q = wm.var[:nw][n][:q][a]
    h_i = wm.var[:nw][n][:h][i]
    h_j = wm.var[:nw][n][:h][j]

    # Return the variables.
    return q, h_i, h_j
end

"Get variables and constants used in the construction of Darcy-Weisbach constraints."
function get_dw_requirements{T <: AbstractWaterFormulation}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    q, h_i, h_j = get_common_variables(wm, a, n)
    viscosity = wm.ref[:nw][n][:options]["viscosity"]
    lambda = calc_friction_factor_dw(wm.ref[:nw][n][:pipes][a], viscosity)
    return q, h_i, h_j, viscosity, lambda
end

"Get variables and constants used in the construction of Hazen-Williams constraints."
function get_hw_requirements{T <: AbstractWaterFormulation}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    q, h_i, h_j = get_common_variables(wm, a, n)
    lambda = calc_friction_factor_hw(wm.ref[:nw][n][:pipes][a])
    return q, h_i, h_j, lambda
end

"These problem forms use binary variables to specify flow direction."
AbstractRelaxedForm = Union{AbstractMICPForm, AbstractMILPRForm}

"Create variables associated with the head for the MICP and MILP-R problems."
function variable_head{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    # Create head variables.
    variable_head_common(wm, n)

    # Create variables that correspond to the absolute value of the head difference.
    diff_min, diff_max = calc_head_difference_bounds(wm.ref[:nw][n][:pipes])
    wm.var[:nw][n][:gamma] = @variable(wm.model, [id in keys(wm.ref[:nw][n][:pipes])],
                                       lowerbound = diff_min[id], upperbound = diff_max[id],
                                       basename = "gamma_$(n)")

    # Create variables that correspond to flow moving from i to j.
    wm.var[:nw][n][:yp] = @variable(wm.model, [id in keys(wm.ref[:nw][n][:pipes])],
                                    category = :Bin, basename = "yp_$(n)")

    # Create variables that correspond to flow moving from j to i.
    wm.var[:nw][n][:yn] = @variable(wm.model, [id in keys(wm.ref[:nw][n][:pipes])],
                                    category = :Bin, basename = "yn_$(n)")
end

"Create variables associated with the head for the MICP and MILP-R problems."
function variable_head_ne{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    # Create head variables.
    variable_head_common(wm, n)

    # Create variables that correspond to flow moving from i to j.
    wm.var[:nw][n][:yp] = @variable(wm.model, [id in keys(wm.ref[:nw][n][:pipes])],
                                    category = :Bin, basename = "yp_$(n)")

    # Create variables that correspond to flow moving from j to i.
    wm.var[:nw][n][:yn] = @variable(wm.model, [id in keys(wm.ref[:nw][n][:pipes])],
                                    category = :Bin, basename = "yn_$(n)")
end

"Create variables associated with the pipe for the MICP and MILP-R problems."
function variable_pipe_ne{T <: AbstractRelaxedForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    # Create pipe variables.
    variable_pipe_ne_common(wm, n)

    # Set up required data to initialize junction variables.
    wm.var[:nw][n][:wp] = Dict{String, Any}()
    wm.var[:nw][n][:wn] = Dict{String, Any}()

    for (pipe_id, pipe) in wm.ref[:nw][n][:ne_pipe]
        # Create binary variables associated with the combination of whether a
        # diameter is used and the direction of flow in a pipe.
        diameters = [d["diameter"] for d in pipe["diameters"]]
        wm.var[:nw][n][:wp][pipe_id] = @variable(wm.model, [d in diameters], category = :Bin,
                                                basename = "wp_$(n)_$(pipe_id)", start = 0)
        wm.var[:nw][n][:wn][pipe_id] = @variable(wm.model, [d in diameters], category = :Bin,
                                                basename = "wn_$(n)_$(pipe_id)", start = 0)
    end
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
    gamma_lb = getlowerbound(h_i) - getupperbound(h_j)
    gamma_ub = getupperbound(h_i) - getlowerbound(h_j)

    # Add the required constraints to define gamma.
    @constraint(wm.model, h_j - h_i + gamma_lb * (y_p - y_n + 1) <= gamma)
    @constraint(wm.model, h_i - h_j + gamma_ub * (y_p - y_n - 1) <= gamma)
    @constraint(wm.model, h_j - h_i + gamma_ub * (y_p - y_n + 1) >= gamma)
    @constraint(wm.model, h_i - h_j + gamma_lb * (y_p - y_n - 1) >= gamma)

    # Get the sum of all junction demands.
    junctions = values(wm.ref[:nw][n][:junctions])
    sum_demand = sum(junction["demand"] for junction in junctions)

    # Add the required constraints to define y_p and y_n.
    @constraint(wm.model, (y_p - 1) * sum_demand <= q)
    @constraint(wm.model, (1 - y_n) * sum_demand >= q)
    @constraint(wm.model, (1 - y_p) * gamma_lb <= h_i - h_j)
    @constraint(wm.model, (1 - y_n) * gamma_ub >= h_i - h_j)
    @constraint(wm.model, y_p + y_n == 1)
end

"These problem forms use exact versions of the head loss equation when the flow direction is fixed."
AbstractUndirectedForm = Union{AbstractMILPForm, AbstractMINLPBForm, AbstractNLPForm}

"Create variables associated with the head for forms of the problem without direction variables."
function variable_head{T <: AbstractUndirectedForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    variable_head_common(wm, n)
end

"Create variables associated with the head for forms of the problem without direction variables."
function variable_head_ne{T <: AbstractUndirectedForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    variable_head_common(wm, n)
end

"These problem forms use exact versions of the head loss equation when the flow direction is fixed."
AbstractExactForm = Union{AbstractMINLPBForm, AbstractNLPForm}
AbstractEqualityForm = Union{AbstractMILPForm, AbstractMINLPBForm, AbstractNLPForm}

"Create variables associated with the pipe for the MICP and MILP-R problems."
function variable_pipe_ne{T <: AbstractEqualityForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    # Create pipe variables.
    variable_pipe_ne_common(wm, n)
end

"Exact Hazen-Williams constraint for flow with known direction."
function constraint_hw_known_direction{T <: AbstractExactForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a non-convex constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) == lambda * (q^2)^0.926)
end

"Exact Darcy-Weisbach constraint with known direction."
function constraint_dw_known_direction{T <: AbstractExactForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a non-convex quadratic constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) == lambda * q^2)
end

function set_initial_solution_ne{T <: AbstractEqualityForm}(wm::GenericWaterModel{T}, wm_solved::GenericWaterModel)
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

        diameter_sol = wm_solved.ref[:nw][wm_solved.cnw][:pipes][ij]["diameter"]
        diameters = [key[1] for key in keys(wm.var[:nw][wm.cnw][:psi][ij])]

        for diameter in diameters
            if diameter == diameter_sol
                pipe = wm.ref[:nw][wm.cnw][:ne_pipe][ij]
                pipe_sol = wm_solved.ref[:nw][wm_solved.cnw][:pipes][ij]
                lambda = calc_friction_factor_hw_ne(pipe_sol, diameter)
                gamma = (h_i_sol - h_j_sol) / lambda
                cost_per_unit_length = [d["costPerUnitLength"] for d in filter((d) -> d["diameter"] == diameter_sol, pipe["diameters"])][1]
                objective_value += pipe["length"] * cost_per_unit_length
                setvalue(wm.var[:nw][wm.cnw][:gamma][ij][diameter], gamma)
                setvalue(wm.var[:nw][wm.cnw][:gamma_sum][ij], gamma)
                setvalue(wm.var[:nw][wm.cnw][:psi][ij][diameter], 1)
            else
                setvalue(wm.var[:nw][wm.cnw][:gamma][ij][diameter], 0.0)
                setvalue(wm.var[:nw][wm.cnw][:psi][ij][diameter], 0)
            end
        end
    end

    setvalue(wm.var[:nw][wm.cnw][:objective], objective_value * 1.0e-6)
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
        setvalue(wm.var[:nw][wm.cnw][:yp][ij], 1 * (Int(flow_direction) >= 0))
        setvalue(wm.var[:nw][wm.cnw][:yn][ij], 1 * (Int(flow_direction) < 0))

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


function constraint_junction_mass_flow{T <: AbstractEqualityForm}(wm::GenericWaterModel{T}, i, n::Int = wm.cnw)
    constraint_flow_conservation(wm, i, n)
end

"Constraints used to define the head difference in the MILP, MINLP-B, and NLP expansion planning problems."
function constraint_define_gamma_hw_ne{T <: AbstractEqualityForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

    # Add a constraint that says at most one diameter must be selected.
    @constraint(wm.model, sum(wm.var[:nw][n][:psi][a]) == 1)

    # Get the pipe associated with the pipe index a.
    pipe = wm.ref[:nw][n][:ne_pipe][a]

    # Get the various pipe diameters from which we can select.
    diameters = [key[1] for key in keys(wm.var[:nw][n][:psi][a])]

    # Constrain each gamma variable.
    for diameter in diameters
        # Gather the required variables and constants.
        psi_d = wm.var[:nw][n][:psi][a][diameter]
        gamma_d = wm.var[:nw][n][:gamma][a][diameter]
        lambda_d = calc_friction_factor_hw_ne(pipe, diameter)
        gamma_ub = (getupperbound(h_i) - getlowerbound(h_j)) / lambda_d
        gamma_lb = (getlowerbound(h_i) - getupperbound(h_j)) / lambda_d

        # Add the four required McCormick constraints.
        @constraint(wm.model, psi_d * gamma_lb <= gamma_d)
        @constraint(wm.model, psi_d * gamma_ub >= gamma_d)
        @constraint(wm.model, (h_i - h_j) / lambda_d - (1 - psi_d) * gamma_ub <= gamma_d)
        @constraint(wm.model, (h_i - h_j) / lambda_d - (1 - psi_d) * gamma_lb >= gamma_d)
    end
end

"These problem forms use different variables for head loss when diameters can be varied."
AbstractInequalityForm = Union{AbstractMICPForm, AbstractMILPRForm}

"Constraints used to define the head difference in the MICP and MILP-R expansion planning problems."
function constraint_define_gamma_hw_ne{T <: AbstractInequalityForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)
    y_p = wm.var[:nw][n][:yp][a]
    y_n = wm.var[:nw][n][:yn][a]

    # Add a constraint that says at most one diameter must be selected.
    @constraint(wm.model, sum(wm.var[:nw][n][:psi][a]) == 1)

    # Get the sum of all junction demands.
    junctions = values(wm.ref[:nw][n][:junctions])
    sum_demand = sum(junction["demand"] for junction in junctions)
    gamma_sum_lb = getlowerbound(wm.var[:nw][n][:gamma_sum][a])
    gamma_sum_ub = getupperbound(wm.var[:nw][n][:gamma_sum][a])
    @constraint(wm.model, sum(wm.var[:nw][n][:wp][a]) == y_p)
    @constraint(wm.model, sum(wm.var[:nw][n][:wn][a]) == y_n)
    @constraint(wm.model, y_p + y_n == 1)
    @constraint(wm.model, (y_p - 1) * sum_demand <= q)
    @constraint(wm.model, (1 - y_n) * sum_demand >= q)
    #@constraint(wm.model, (1 - y_p) * gamma_sum_lb <= h_i - h_j)
    #@constraint(wm.model, (1 - y_n) * gamma_sum_ub >= h_i - h_j)

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

        # Create an auxiliary variable for the product of psi_d and (y_p - y_n).
        #@constraint(wm.model, (z_p - z_n) + psi_d >= 0)
        #@constraint(wm.model, -(z_p - z_n) + psi_d >= 0)
        #@constraint(wm.model, -(z_p - z_n) - psi_d + (y_p - y_n) + 1 >= 0)
        #@constraint(wm.model, (z_p - z_n) - psi_d - (y_p - y_n) + 1 >= 0)

        #@constraint(wm.model, (z_p - 1) * sum_demand <= q)
        #@constraint(wm.model, (1 - z_n) * sum_demand >= q)
        #@constraint(wm.model, (1 - z_p) * gamma_lb <= h_i - h_j)
        #@constraint(wm.model, (1 - z_n) * gamma_ub >= h_i - h_j)

        # Ensure that z is chosen when the diameter is chosen.
        @constraint(wm.model, (z_p + z_n) == psi_d)

        # Link z_p and z_n with q and (h_i - h_j).
        @constraint(wm.model, (z_p - 1) * sum_demand <= q)
        @constraint(wm.model, (1 - z_n) * sum_demand >= q)
        @constraint(wm.model, (1 - z_p) * gamma_lb <= (h_i - h_j) / lambda_d)
        @constraint(wm.model, (1 - z_n) * gamma_ub >= (h_i - h_j) / lambda_d)

        # McCormick relaxation of (z_p - z_n) * (h_i - h_j) / lambda_d == gamma_d.
        #@NLconstraint(wm.model, (z_p - z_n) * (h_i - h_j) / lambda_d == gamma_d)

        h_diff = @variable(wm.model, lowerbound = gamma_lb, upperbound = gamma_ub) 
        z_diff = @variable(wm.model, category = :Int, lowerbound = -1, upperbound = 1)
        @constraint(wm.model, h_diff == (h_i - h_j) / lambda_d)
        @constraint(wm.model, z_diff == z_p - z_n)
        InfrastructureModels.relaxation_product(wm.model, h_diff, z_diff, gamma_d)

        ## Create McCormick constraints that ensure gamma_d is positive.
        #@constraint(wm.model, (h_j - h_i) / lambda_d + gamma_lb * (z_p - z_n + 1) <= gamma_d)
        #@constraint(wm.model, (h_i - h_j) / lambda_d + gamma_ub * (z_p - z_n - 1) <= gamma_d)
        #@constraint(wm.model, (h_j - h_i) / lambda_d + gamma_ub * (z_p - z_n + 1) >= gamma_d)
        #@constraint(wm.model, (h_i - h_j) / lambda_d + gamma_lb * (z_p - z_n - 1) >= gamma_d)
    end
end

"Constraint to ensure at least one direction is set to take flow away from a source."
function constraint_source_flow{T <: AbstractInequalityForm}(wm::GenericWaterModel{T}, i, f_branches, t_branches, n::Int = wm.cnw)
    y_p = wm.var[:nw][n][:yp]
    y_n = wm.var[:nw][n][:yn]
    @constraint(wm.model, sum(y_p[a] for a in f_branches) + sum(y_n[a] for a in t_branches) >= 1)
end

"Constraint to ensure at least one direction is set to take flow to a junction with demand."
function constraint_sink_flow{T <: AbstractInequalityForm}(wm::GenericWaterModel{T}, i, f_branches, t_branches, n::Int = wm.cnw)
    y_p = wm.var[:nw][n][:yp]
    y_n = wm.var[:nw][n][:yn]
    @constraint(wm.model, sum(y_n[a] for a in f_branches) + sum(y_p[a] for a in t_branches) >= 1)
end

function constraint_junction_mass_flow{T <: AbstractInequalityForm}(wm::GenericWaterModel{T}, i, n::Int = wm.cnw)
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

function constraint_no_good_ne{T <: AbstractInequalityForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    yp_solution = getvalue(wm.var[:nw][n][:yp])
    yp_ones = Array{JuMP.Variable}([wm.var[:nw][n][:yp][idx[1]] for idx in keys(yp_solution) if yp_solution[idx[1]] > 1.0e-4])
    yp_zeros = Array{JuMP.Variable}([wm.var[:nw][n][:yp][idx[1]] for idx in keys(yp_solution) if yp_solution[idx[1]] <= 1.0e-4])

    pipe_ids = ids(wm, :ne_pipe)
    psi_vars = Array{JuMP.Variable}([wm.var[:nw][n][:psi][a][d["diameter"]] for a in ids(wm, :ne_pipe) for d in wm.ref[:nw][n][:ne_pipe][a]["diameters"]])
    psi_ones = Array{JuMP.Variable}([psi for psi in psi_vars if getvalue(psi) > 1.0e-4])
    psi_zeros = Array{JuMP.Variable}([psi for psi in psi_vars if getvalue(psi) <= 1.0e-4])

    one_vars = psi_ones #vcat(yp_ones, psi_ones)
    zero_vars = psi_zeros #vcat(yp_zeros, psi_zeros)

    @constraint(wm.model, sum(zero_vars) - sum(one_vars) >= 1 - length(one_vars))
end

function constraint_degree_two{T <: AbstractInequalityForm}(wm::GenericWaterModel{T}, idx, n::Int = wm.cnw)
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
