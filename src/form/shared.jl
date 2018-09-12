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
    @constraint(wm.model, (1 - y_p) * (h_i_lb - h_j_ub) <= h_i - h_j)
    @constraint(wm.model, (1 - y_n) * (h_i_ub - h_j_lb) >= h_i - h_j)
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

"These problem forms use different variables for head loss when diameters can be varied."
AbstractEqualityForm = Union{AbstractMILPForm, AbstractMINLPBForm, AbstractNLPForm}

"Constraints used to define the head difference in the MILP, MINLP-B, and NLP expansion planning problems."
function constraint_define_gamma_hw_ne{T <: AbstractEqualityForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j = get_common_variables(wm, a, n)

    # Add a constraint that says at most one diameter must be selected.
    @constraint(wm.model, sum(wm.var[:nw][n][:psi][a]) <= 1)

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
    @constraint(wm.model, y_p + y_n == 1)

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

        # Create an auxiliary variable for the product of psi_d and (y_p - y_n).
        w_d = @variable(wm.model, category = :Int, lowerbound = -1, upperbound = 1, start = 0)
        @constraint(wm.model, w_d + psi_d >= 0)
        @constraint(wm.model, -w_d + psi_d >= 0)
        @constraint(wm.model, -w_d - psi_d + (y_p - y_n) + 1 >= 0)
        @constraint(wm.model, w_d - psi_d - (y_p - y_n) + 1 >= 0)

        # Add the four required McCormick constraints.
        @constraint(wm.model, gamma_d - gamma_lb * w_d >= 0)
        @constraint(wm.model, -gamma_d + gamma_ub * w_d >= 0)
        @constraint(wm.model, -gamma_d + gamma_lb * w_d + (h_i - h_j) / lambda_d - gamma_lb >= 0)
        @constraint(wm.model, gamma_d - gamma_ub * w_d - (h_i - h_j) / lambda_d + gamma_ub >= 0)
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

    #elseif !has_supply && !has_deman

    #if !has_supply
    #        
    #if fgfirm == 0.0 && flfirm == 0.0 && junction["degree"] == 2
    #    constraint_conserve_flow(wm, i, n)
    #end
end
