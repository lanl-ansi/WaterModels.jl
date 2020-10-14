# Define common CRD (continuous or convex relaxation- and direction-based) implementations
# of water distribution network constraints, which use directed flow variables.


########################################## PIPES ##########################################


function constraint_pipe_head_loss(wm::CRDWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64, L::Float64, r::Float64)
    # Gather directed flow and head difference variables.
    qp, qn = var(wm, n, :qp_pipe, a), var(wm, n, :qn_pipe, a)
    dhp, dhn = var(wm, n, :dhp_pipe, a), var(wm, n, :dhn_pipe, a)

    # Add relaxed constraints for head loss in the positive and negative directions.
    c_1 = JuMP.@NLconstraint(wm.model, r * head_loss(qp) <= inv(L) * dhp)
    c_2 = JuMP.@NLconstraint(wm.model, r * head_loss(qn) <= inv(L) * dhn)

    # Add linear upper bounds on the above convex relaxations.
    rhs_p = r * JuMP.upper_bound(qp)^(exponent - 1.0) * qp
    c_3 = JuMP.@constraint(wm.model, inv(L) * dhp <= rhs_p)
    rhs_n = r * JuMP.upper_bound(qn)^(exponent - 1.0) * qn
    c_4 = JuMP.@constraint(wm.model, inv(L) * dhn <= rhs_n)

    # Append the :pipe_head_loss constraint array.
    append!(con(wm, n, :pipe_head_loss)[a], [c_1, c_2, c_3, c_4])
end


function constraint_on_off_pipe_head_loss_des(wm::CRDWaterModel, n::Int, a::Int, exponent::Float64, node_fr::Int, node_to::Int, L::Float64, resistances)
    # Collect head difference and flow variables.
    dhp, dhn = var(wm, n, :dhp_des_pipe, a), var(wm, n, :dhn_des_pipe, a)
    qp, qn = var(wm, n, :qp_des_pipe, a), var(wm, n, :qn_des_pipe, a)

    for (r_id, r) in enumerate(resistances)
        # Build the relaxed head loss constraints.
        c_1 = JuMP.@NLconstraint(wm.model, r * head_loss(qp[r_id]) <= inv(L) * dhp)
        c_2 = JuMP.@NLconstraint(wm.model, r * head_loss(qn[r_id]) <= inv(L) * dhn)

        # Append the :on_off_pipe_head_loss_des constraint array.
        append!(con(wm, n, :on_off_pipe_head_loss_des, a), [c_1, c_2])
    end

    # Add linear upper bounds for the positive portion of head loss.
    qp_ub = JuMP.upper_bound.(qp)
    slopes_p = resistances .* qp_ub.^(exponent - 1.0)
    c_3 = JuMP.@constraint(wm.model, inv(L)*dhp <= sum(slopes_p .* qp))

    # Add linear upper bounds for the negative portion of head loss.
    qn_ub = JuMP.upper_bound.(qn)
    slopes_n = resistances .* qn_ub.^(exponent - 1.0)
    c_4 = JuMP.@constraint(wm.model, inv(L)*dhn <= sum(slopes_n .* qn))

    # Append the :on_off_pipe_head_loss_des constraint array.
    append!(con(wm, n, :on_off_pipe_head_loss_des)[a], [c_3, c_4])
end


########################################## PUMPS ##########################################


"Instantiate variables associated with modeling each pump's head gain."
function variable_pump_head_gain(wm::CRDWaterModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    # Initialize variables for total hydraulic head gain from a pump.
    g = var(wm, nw)[:g_pump] = JuMP.@variable(wm.model, [a in ids(wm, nw, :pump)],
        base_name="$(nw)_g_pump", lower_bound=0.0, # Pump gain is nonnegative.
        start=comp_start_value(ref(wm, nw, :pump, a), "g_pump_start"))

    # Initialize an entry to the solution component dictionary for head gains.
    report && sol_component_value(wm, nw, :pump, :g, ids(wm, nw, :pump), g)

    # If the number of breakpoints is not positive, return.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 2)

    # Create weights involved in convex combination constraints.
    lambda = var(wm, nw)[:lambda_pump] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pump), k in 1:pump_breakpoints],
        base_name="$(nw)_lambda", lower_bound=0.0, upper_bound=1.0,
        start=comp_start_value(ref(wm, nw, :pump, a), "lambda_start", k))

    # Create binary variables involved in convex combination constraints.
    x_pw = var(wm, nw)[:x_pw_pump] = JuMP.@variable(wm.model,
        [a in ids(wm, nw, :pump), k in 1:pump_breakpoints-1],
        base_name="$(nw)_x_pw", binary=true,
        start=comp_start_value(ref(wm, nw, :pump, a), "x_pw_start", k))
end


"Add constraints associated with modeling a pump's head gain."
function constraint_on_off_pump_head_gain(wm::CRDWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, pc::Array{Float64}, q_min_forward::Float64)
    # Get the number of breakpoints for the pump.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 2)

    # Gather pump flow, head gain, and status variables.
    qp, g, z = var(wm, n, :qp_pump, a), var(wm, n, :g_pump, a), var(wm, n, :z_pump, a)

    # Define the (relaxed) head gain relationship for the pump.
    c_1 = JuMP.@constraint(wm.model, g <= pc[1]*qp^2 + pc[2]*qp + pc[3]*z)

    # Gather flow, head gain, and convex combination variables.
    lambda, x_pw = var(wm, n, :lambda_pump), var(wm, n, :x_pw_pump)

    # Add the required SOS constraints.
    c_2 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == z)
    c_3 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == z)
    c_4 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
    c_5 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

    # Add a constraint for the flow piecewise approximation.
    breakpoints = range(q_min_forward, stop=JuMP.upper_bound(qp), length=pump_breakpoints)
    qp_lhs = sum(breakpoints[k] * lambda[a, k] for k in 1:pump_breakpoints)
    c_6 = JuMP.@constraint(wm.model, qp_lhs == qp)

    # Add a constraint that lower-bounds the head gain variable.
    f = (pc[1] .* breakpoints.^2) .+ (pc[2] .* breakpoints) .+ pc[3]
    gain_lb_expr = sum(f[k] .* lambda[a, k] for k in 1:pump_breakpoints)
    c_7 = JuMP.@constraint(wm.model, gain_lb_expr <= g)

    # Append the constraint array.
    append!(con(wm, n, :on_off_pump_head_gain, a), [c_1, c_2, c_3, c_4, c_5, c_6, c_7])

    for k in 2:pump_breakpoints-1
        # Add the adjacency constraints for piecewise variables.
        adjacency = x_pw[a, k-1] + x_pw[a, k]
        c_8_k = JuMP.@constraint(wm.model, lambda[a, k] <= adjacency)
        append!(con(wm, n, :on_off_pump_head_gain, a), [c_8_k])
    end
end


######################################## OBJECTIVES ########################################


"Instantiate the objective associated with the Optimal Water Flow problem."
function objective_owf(wm::CRDWaterModel)
    # Get the number of breakpoints for the pump.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 2)

    # Initialize the objective function.
    objective = JuMP.AffExpr(0.0)

    for (n, nw_ref) in nws(wm)
        # Get common variables.
        lambda = var(wm, n, :lambda_pump)

        # Get common constant parameters.
        constant = _DENSITY * _GRAVITY * ref(wm, n, :time_step)

        for (a, pump) in nw_ref[:pump]
            if haskey(pump, "energy_price")
                # Get price and pump curve data.
                price = pump["energy_price"]
                head_curve = ref(wm, n, :pump, a)["head_curve"]
                curve_fun = _get_function_from_head_curve(head_curve)

                # Get flow-related variables and data.
                qp, z = var(wm, n, :qp_pump, a), var(wm, n, :z_pump, a)
                qp_ub = JuMP.upper_bound(qp)

                # Generate a set of uniform flow and cubic function breakpoints.
                breakpoints = range(0.0, stop=qp_ub, length=pump_breakpoints)
                f = _calc_cubic_flow_values(collect(breakpoints), curve_fun)

                # Get pump efficiency data.
                if haskey(pump, "efficiency_curve")
                    eff_curve = pump["efficiency_curve"]
                    eff = _calc_efficiencies(collect(breakpoints), eff_curve)
                else
                    eff = pump["efficiency"]
                end

                # Add the cost corresponding to the current pump's operation.
                inner_expr = (constant*price) .* inv.(eff) .* f
                cost = sum(inner_expr[k]*lambda[a, k] for k in 1:pump_breakpoints)
                JuMP.add_to_expression!(objective, cost)
            else
                Memento.error(_LOGGER, "No cost given for pump \"$(pump["name"])\"")
            end
        end
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end
