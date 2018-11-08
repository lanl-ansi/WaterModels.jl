"These problem forms use exact versions of the head loss equation when the flow direction is fixed."
AbstractEqualityForm = Union{AbstractMILPForm, AbstractMINLPBForm, AbstractNLPForm}

"Create variables associated with the head for forms of the problem without direction variables."
function variable_head{T <: AbstractEqualityForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    variable_head_common(wm, n)
end

"Create variables associated with the head for forms of the problem without direction variables."
function variable_head_ne{T <: AbstractEqualityForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    variable_head_common(wm, n)
end

"Create variables associated with the pipe for the MICP and MILP-R problems."
function variable_pipe_ne{T <: AbstractEqualityForm}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    variable_pipe_ne_common(wm, n)
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
