# Define MILPR (mixed-integer linear, relaxed program) implementations of water distribution models.

function get_objective_values(breakpoints::Array{Float64}, curve_fun::Array{Float64})
    return [curve_fun[1]*x^3 + curve_fun[2]*x^2 + curve_fun[3]*x for x in breakpoints]
end

function get_cubic_outer_approximation(q::JuMP.VariableRef, x_pump::JuMP.VariableRef, q_hat::Float64, curve_fun::Array{Float64})
    return curve_fun[1]*q_hat^3*x_pump + 3.0*curve_fun[1] * q_hat^2 * (q - q_hat)*x_pump
end

function get_obj_linear_outer_approximation(q::JuMP.VariableRef, x_pump::JuMP.VariableRef, q_hat::Float64, curve_fun::Array{Float64})
    f_a = curve_fun[1]*q_hat^3 + curve_fun[2]*q_hat^2 + curve_fun[3]*q_hat
    return f_a + 3.0*curve_fun[1] * (q - q_hat)^2 + 2.0*curve_fun[2] * (q - q_hat) + curve_fun[3]
end

function get_linear_outer_approximation(q::JuMP.VariableRef, q_hat::Float64, alpha::Float64)
    return q_hat^alpha + alpha * q_hat^(alpha - 1.0) * (q - q_hat)
end

function get_linear_outer_approximation_pump(q::JuMP.VariableRef, x_pump::JuMP.VariableRef, q_hat::Float64, a::Float64, b::Float64, c::Float64)
    return a*q_hat^2 + b*q_hat + c + 2.0*a*(q - q_hat) + b
end

function get_linear_outer_approximation_pump_on(q::JuMP.VariableRef, q_hat::Float64, a::Float64, b::Float64, c::Float64)
    return a*q_hat^2 + b*q_hat + c + 2.0*a*(q - q_hat) + b
end

function variable_flow_piecewise_weights(wm::AbstractMILPRModel, n::Int=wm.cnw)
    # Set the number of breakpoints used in each outer-approximation.
    num_breakpoints = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 2

    var(wm, n)[:lambda] = JuMP.@variable(wm.model, [a in ids(wm, n, :pump),
        k in 1:num_breakpoints], base_name="lambda[$(n)]", lower_bound=0.0,
        upper_bound=1.0, start=get_start(ref(wm, n, :link), a, "lambda_start", 0.0))
end

function variable_flow_piecewise_adjacency(wm::AbstractMILPRModel, n::Int=wm.cnw)
    # Set the number of breakpoints used in each outer-approximation.
    num_breakpoints = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 2
    var(wm, n)[:x_pw] = JuMP.@variable(wm.model, [a in ids(wm, n, :pump),
        k in 1:num_breakpoints-1], base_name="x_pw[$(n)]", binary=true,
        start=get_start(ref(wm, n, :link_fixed), a, "x_pw_start", 0.0))
end

function variable_flow(wm::AbstractMILPRModel; nw::Int=wm.cnw, bounded::Bool=true)
    # Create directed flow variables (i.e., qp and qn).
    variable_flow_common(wm, nw=nw, bounded=true)

    # Variables for direction.
    variable_flow_direction(wm, nw=nw)

    ## Create variables required for objective's piecewise approximation.
    #variable_flow_piecewise_weights(wm, nw)
    #variable_flow_piecewise_adjacency(wm, nw)
end

"Pump head gain constraint when the pump status is ambiguous."
function constraint_head_gain_pump(wm::AbstractMILPRModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Gather common variables.
    q = var(wm, n, :qp, a)
    dhp = var(wm, n, :dhp, a)
    dhn = var(wm, n, :dhn, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)
    x_pump = var(wm, n, :x_pump, a)
    q_ub = JuMP.has_upper_bound(q) ? JuMP.upper_bound(q) : 10.0

    # Get head difference lower bounds.
    dh_lb = JuMP.lower_bound(h_j) - JuMP.upper_bound(h_i)
    dh_ub = JuMP.upper_bound(h_j) - JuMP.lower_bound(h_i)

    # If the pump is off, decouple the head difference relationship.
    g = curve_fun[1]*q*q + curve_fun[2]*q + curve_fun[3]*x_pump
    c_1 = JuMP.@constraint(wm.model, (h_j - h_i) - g >= dh_lb * (1.0 - x_pump))
    c_2 = JuMP.@constraint(wm.model, (h_j - h_i) - g <= dh_ub * (1.0 - x_pump))

    # If the pump is off, the flow along the pump must be zero.
    c_3 = JuMP.@constraint(wm.model, q <= q_ub * x_pump)
    c_4 = JuMP.@constraint(wm.model, q >= 6.31465679e-6 * x_pump)

    # Append the constraint array.
    append!(con(wm, n, :head_gain)[a], [c_1, c_2, c_3, c_4])
end

"Pump head gain constraint when the pump is forced to be on."
function constraint_head_gain_pump_on(wm::AbstractMILPRModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Fix reverse flow variable to zero (since this is a pump).
    JuMP.fix(var(wm, n, :qn, a), 0.0, force=true)

    # Gather common variables.
    qp = var(wm, n, :qp, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)

    # Define the (relaxed) head gain caused by the pump.
    qp_ub = JuMP.has_upper_bound(qp) ? JuMP.upper_bound(qp) : 10.0
    num_breakpoints = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 1
    breakpoints = range(0.0, stop=qp_ub, length=num_breakpoints)

    for q_hat in breakpoints
        cut_lhs = get_linear_outer_approximation_pump_on(qp, q_hat, curve_fun[1], curve_fun[2], curve_fun[3])
        c = JuMP.@constraint(wm.model, cut_lhs <= h_j - h_i)
        append!(con(wm, n, :head_gain)[a], [c])
    end
end

function constraint_head_loss_pipe_ne(wm::AbstractMILPRModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, resistances)
    # Set the number of breakpoints used in each outer-approximation.
    num_breakpoints = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 1

    for (r_id, r) in enumerate(resistances)
        qp_ne = var(wm, n, :qp_ne, a)[r_id]
        qp_ne_ub = JuMP.has_upper_bound(qp_ne) ? JuMP.upper_bound(qp_ne) : 10.0
        dhp = var(wm, n, :dhp, a)

        if qp_ne_ub > 0.0 && num_breakpoints > 0
            breakpoints = range(0.0, stop=qp_ne_ub, length=num_breakpoints+2)

            for q_hat in breakpoints[2:num_breakpoints+1]
                cut_lhs = r * get_linear_outer_approximation(qp_ne, q_hat, alpha)
                con_p = JuMP.@constraint(wm.model, cut_lhs <= inv(L) * dhp)
            end
        end

        qn_ne = var(wm, n, :qn_ne, a)[r_id]
        qn_ne_ub = JuMP.has_upper_bound(qn_ne) ? JuMP.upper_bound(qn_ne) : 10.0
        dhn = var(wm, n, :dhn, a)

        if qn_ne_ub > 0.0 && num_breakpoints > 0
            breakpoints = range(0.0, stop=qn_ne_ub, length=num_breakpoints+2)

            for q_hat in breakpoints[2:num_breakpoints+1]
                cut_lhs = r * get_linear_outer_approximation(qn_ne, q_hat, alpha)
                con_n = JuMP.@constraint(wm.model, cut_lhs <= inv(L) * dhn)
            end
        end
    end
end

function constraint_head_loss_pipe(wm::AbstractMILPRModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    # Set the number of breakpoints used in each outer-approximation.
    num_breakpoints = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 1

    qp = var(wm, n, :qp, a)
    qp_ub = JuMP.has_upper_bound(qp) ? JuMP.upper_bound(qp) : 10.0
    dhp = var(wm, n, :dhp, a)

    if qp_ub > 0.0 && num_breakpoints > 0
        for q_hat in range(0.0, stop=qp_ub, length=num_breakpoints)
            cut_lhs = r * get_linear_outer_approximation(qp, q_hat, alpha)
            con_p = JuMP.@constraint(wm.model, cut_lhs <= inv(L) * dhp)
        end
    end

    qn = var(wm, n, :qn, a)
    qn_ub = JuMP.has_upper_bound(qn) ? JuMP.upper_bound(qn) : 10.0
    dhn = var(wm, n, :dhn, a)

    if qn_ub > 0.0 && num_breakpoints > 0
        for q_hat in range(0.0, stop=qn_ub, length=num_breakpoints)
            cut_lhs = r * get_linear_outer_approximation(qn, q_hat, alpha)
            con_n = JuMP.@constraint(wm.model, cut_lhs <= inv(L) * dhn)
        end
    end
end

function constraint_head_loss_check_valve(wm::AbstractMILPRModel, n::Int, a::Int, node_fr::Int, node_to::Int, L::Float64, r::Float64) 
    # Set the number of breakpoints used in each outer-approximation.
    num_breakpoints = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 1
    alpha = ref(wm, n, :alpha)
    x_cv = var(wm, n, :x_cv, a)

    qp = var(wm, n, :qp, a)
    qp_ub = JuMP.has_upper_bound(qp) ? JuMP.upper_bound(qp) : 10.0
    dhp = var(wm, n, :dhp, a)

    if qp_ub > 0.0 && num_breakpoints > 0
        for q_hat in range(0.0, stop=qp_ub, length=num_breakpoints)
            cut_lhs = get_linear_outer_approximation(qp, q_hat, alpha)
            con_p = JuMP.@constraint(wm.model, r * cut_lhs <= inv(L) * dhp)
            append!(con(wm, n, :head_loss)[a], [con_p])
        end
    end

    dhn = var(wm, n, :dhn, a)
    dhn_ub = JuMP.upper_bound(dhn)
    con_n = JuMP.@constraint(wm.model, dhn <= dhn_ub * (1.0 - x_cv))
end

function objective_wf(wm::AbstractMILPRModel)
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end

function objective_owf(wm::AbstractMILPRModel)
    objective = JuMP.QuadExpr()
    num_breakpoints = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 2

    for (n, nw_ref) in nws(wm)
        efficiency = 0.85 # TODO: Change this after discussion. 0.85 follows Fooladivanda.
        rho = 1000.0 # Water density (kilogram per cubic meter).
        g = 9.81 # Gravitational acceleration (meter per second squared).
        time_step = nw_ref[:option]["time"]["hydraulic_timestep"]
        constant = g*rho*time_step*inv(efficiency)

        for (a, pump) in nw_ref[:pump]
            cubic_var = JuMP.@variable(wm.model)
            qp = var(wm, n, :qp, a)
            #lambda = var(wm, n, :lambda)
            #x_pw = var(wm, n, :x_pw)
            x_pump = var(wm, n, :x_pump, a)

            #c_1 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == x_pump)
            #c_2 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == x_pump)
            #c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
            #c_4 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

            pump_curve = ref(wm, n, :pump, a)["pump_curve"]
            curve_fun = _get_function_from_pump_curve(pump_curve)
            qp_ub = JuMP.has_upper_bound(qp) ? JuMP.upper_bound(qp) : 10.0
            breakpoints = range(0.0, stop=qp_ub, length=num_breakpoints)

            for q_hat in breakpoints
                lhs = get_cubic_outer_approximation(qp, x_pump, q_hat, curve_fun)
                JuMP.@constraint(wm.model, lhs <= cubic_var)
            end

            energy_price = 0.01 * pump["energy_price"]
            quad_terms = curve_fun[2]*qp^2 + curve_fun[3]*qp
            JuMP.add_to_expression!(objective, constant*energy_price * (cubic_var + quad_terms))

            #f = get_objective_values(collect(breakpoints), curve_fun)
            #println(curve_fun)
            #lhs = energy_price*constant * sum(f[k] .* lambda[a, k] for k in 1:num_breakpoints)
            #q_lhs = sum(breakpoints[k] * lambda[a, k] for k in 1:num_breakpoints)
            #c_6 = JuMP.@constraint(wm.model, q_lhs == qp)

            #for k in 2:num_breakpoints-1
            #    c_7 = JuMP.@constraint(wm.model, lambda[a, k] <= x_pw[a, k-1] + x_pw[a, k])
            #end

            #JuMP.add_to_expression!(objective, lhs)
        end
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end

#function objective_owf(wm::AbstractMILPRModel) 
#    objective = JuMP.AffExpr(0.0)
#
#    for (n, nw_ref) in nws(wm)
#        costs = JuMP.@variable(wm.model, [a in ids(wm, n, :pump)],
#            base_name="costs[$(n)]", lower_bound=0.0,
#            start=get_start(ref(wm, n, :pump), a, "costs", 0.0))
#
#        efficiency = 0.85 # TODO: Change this after discussion. 0.85 follows Fooladivanda.
#        rho = 1000.0 # Water density (kilogram per cubic meter).
#        gravity = 9.80665 # Gravitational acceleration (meter per second squared).
#        time_step = 1.0 #nw_ref[:option]["time"]["hydraulic_timestep"]
#        constant = rho * gravity * time_step * inv(efficiency)
#
#        for (a, pump) in nw_ref[:pump]
#            if haskey(pump, "energy_price")
#                energy_price = pump["energy_price"]
#                pump_curve = ref(wm, n, :pump, a)["pump_curve"]
#                curve_fun = _get_function_from_pump_curve(pump_curve)
#
#                qp = var(wm, n)[:qp][a]
#                x_pump = var(wm, n)[:x_pump][a]
#                qp_ub = JuMP.upper_bound(qp)
#
#                num_breakpoints = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 1
#                breakpoints = range(0.0, stop=qp_ub, length=num_breakpoints)
#
#                for qp_hat in breakpoints
#                    rhs = get_obj_linear_outer_approximation(qp, x_pump, qp_hat, curve_fun)
#                    JuMP.@constraint(wm.model, costs[a] >= constant*energy_price * rhs)
#                end
#            else
#                Memento.error(_LOGGER, "No cost given for pump \"$(pump["name"])\"")
#            end
#        end
#
#        JuMP.add_to_expression!(objective, sum(costs))
#    end
#
#    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
#end
