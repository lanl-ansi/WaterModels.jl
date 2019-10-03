# Define MILPR (mixed-integer linear, relaxed program) implementations of water distribution models.

# Define MILPR (non-convex nonlinear programming) implementations of water distribution models.
function variable_pump(wm::AbstractMILPRModel, n::Int=wm.cnw) 
    variable_fixed_speed_pump_operation(wm, n)
    # TODO: The below need not be created in the OWF.
    variable_fixed_speed_pump_threshold(wm, n) 
end

function get_linear_outer_approximation(q::JuMP.VariableRef, q_hat::Float64, alpha::Float64)
    return q_hat^alpha + alpha * q_hat^(alpha - 1.0) * (q - q_hat)
end

"Pump head gain constraint when the pump status is ambiguous."
function constraint_head_gain_pump(wm::AbstractMILPRModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Fix reverse flow variable to zero (since this is a pump).
    qn = var(wm, n, :qn, a)
    JuMP.fix(qn, 0.0, force=true)

    # Gather common variables.
    qp = var(wm, n, :qp, a)
    qp_ub = JuMP.has_upper_bound(qp) ? JuMP.upper_bound(qp) : 10.0

    g = var(wm, n, :g, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)
    x_pump = var(wm, n, :x_pump, a)

    # Define the (relaxed) head gain caused by the pump.
    c_1 = JuMP.@NLconstraint(wm.model, curve_fun[1]*qp^2 + curve_fun[2]*qp + curve_fun[3]*x_pump <= g)

    # If the pump is off, decouple the head difference relationship.
    c_2 = JuMP.@constraint(wm.model, h_j - h_i - g <= 1.0e6 * (1 - x_pump))
    c_3 = JuMP.@constraint(wm.model, h_j - h_i - g >= -1.0e6 * (1 - x_pump))

    # If the pump is off, the flow along the pump must be zero.
    c_4 = JuMP.@constraint(wm.model, qp <= qp_ub * x_pump)
    c_5 = JuMP.@constraint(wm.model, qp >= 1.0e-6 * x_pump)

    # Append the constraint array.
    con(wm, n, :head_gain)[a] = [c_1, c_2, c_3, c_4, c_5]
end

function constraint_head_loss_pipe_ne(wm::AbstractMILPRModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, resistances)
    # Set the number of breakpoints used in each outer-approximation.
    num_breakpoints = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 1

    for (r_id, r) in enumerate(resistances)
        qp_ne = var(wm, n, :qp_ne, a)[r_id]
        qp_ne_ub = JuMP.upper_bound(qp_ne)
        dhp = var(wm, n, :dhp, a)

        if qp_ne_ub > 0.0 && num_breakpoints > 0
            breakpoints = range(0.0, stop=qp_ne_ub, length=num_breakpoints+2)

            for q_hat in breakpoints[2:num_breakpoints+1]
                cut_lhs = r * get_linear_outer_approximation(qp_ne, q_hat, alpha)
                con_p = JuMP.@constraint(wm.model, cut_lhs - inv(L) * dhp <= 0.0)
            end
        end

        qn_ne = var(wm, n, :qn_ne, a)[r_id]
        qn_ne_ub = JuMP.upper_bound(qn_ne)
        dhn = var(wm, n, :dhn, a)

        if qn_ne_ub > 0.0 && num_breakpoints > 0
            breakpoints = range(0.0, stop=qn_ne_ub, length=num_breakpoints+2)

            for q_hat in breakpoints[2:num_breakpoints+1]
                cut_lhs = r * get_linear_outer_approximation(qn_ne, q_hat, alpha)
                con_n = JuMP.@constraint(wm.model, cut_lhs - inv(L) * dhn <= 0.0)
            end
        end
    end
end

function constraint_head_loss_pipe(wm::AbstractMILPRModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    # Set the number of breakpoints used in each outer-approximation.
    num_breakpoints = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 1

    qn = var(wm, n, :qn, a)
    qn_ub = JuMP.upper_bound(qn)
    dhn = var(wm, n, :dhn, a)

    if qn_ub > 0.0 && num_breakpoints > 0
        breakpoints = range(0.0, stop=qn_ub, length=num_breakpoints+2)

        for q_hat in breakpoints[2:num_breakpoints+1]
            cut_lhs = r * get_linear_outer_approximation(qn, q_hat, alpha)
            con_n = JuMP.@constraint(wm.model, cut_lhs - inv(L) * dhn <= 0.0)
        end
    end

    qp = var(wm, n, :qp, a)
    qp_ub = JuMP.upper_bound(qp)
    dhp = var(wm, n, :dhp, a)

    if qp_ub > 0.0 && num_breakpoints > 0
        breakpoints = range(0.0, stop=qp_ub, length=num_breakpoints+2)

        for q_hat in breakpoints[2:num_breakpoints+1]
            cut_lhs = r * get_linear_outer_approximation(qp, q_hat, alpha)
            con_p = JuMP.@constraint(wm.model, cut_lhs - inv(L) * dhp <= 0.0)
        end
    end
end

function objective_wf(wm::AbstractMILPRModel)
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end

function objective_owf(wm::AbstractMILPRModel) 
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end
