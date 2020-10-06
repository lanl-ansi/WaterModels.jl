# Define common CQRD (quadratic relaxation- and direction-based) implementations
# of water distribution network constraints, which use directed flow variables.


function constraint_pipe_head_loss(wm::CQRDWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64, L::Float64, r::Float64)
    # Get the number of breakpoints for the pipe.
    num_breakpoints = get(wm.ext, :pipe_breakpoints, 1)

    # Get the variable for flow directionality.
    y = var(wm, n, :y_pipe, a)

    # Get variables for positive flow and head difference.
    qp, dhp = var(wm, n, :qp_pipe, a), var(wm, n, :dhp_pipe, a)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(0.0, stop = JuMP.upper_bound(qp), length = num_breakpoints+2)[2:end-1]
        # Add a linear outer approximation of the convex relaxation at `pt`.
        lhs = _get_head_loss_oa_binary(qp, y, pt, exponent)
        c = JuMP.@constraint(wm.model, r * lhs <= inv(L) * dhp)
        append!(con(wm, n, :pipe_head_loss)[a], [c])
    end

    # Add linear upper bounds on the above outer approximations.
    rhs = r * JuMP.upper_bound(qp)^(exponent - 1.0) * qp
    c = JuMP.@constraint(wm.model, inv(L) * dhp <= rhs)

    # Append the :pipe_head_loss constraint array.
    append!(con(wm, n, :pipe_head_loss)[a], [c])

    # Get variables for negative flow and head difference.
    qn, dhn = var(wm, n, :qn_pipe, a), var(wm, n, :dhn_pipe, a)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(0.0, stop = JuMP.upper_bound(qn), length = num_breakpoints+2)[2:end-1]
        # Add a linear outer approximation of the convex relaxation at `pt`.
        lhs = _get_head_loss_oa_binary(qn, 1.0 - y, pt, exponent)
        c = JuMP.@constraint(wm.model, r * lhs <= inv(L) * dhn)
        append!(con(wm, n, :pipe_head_loss)[a], [c])
    end

    # Add linear upper bounds on the above outer approximations.
    rhs = r * JuMP.upper_bound(qn)^(exponent - 1.0) * qn
    c = JuMP.@constraint(wm.model, inv(L) * dhn <= rhs)

    # Append the :pipe_head_loss constraint array.
    append!(con(wm, n, :pipe_head_loss)[a], [c])
end


"Pump head gain constraint when the pump status is ambiguous."
function constraint_on_off_pump_head_gain(wm::CQRDWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, pc::Array{Float64}, q_min_active::Float64)
    # Gather pump flow, head gain, and status variables.
    qp, g, z = var(wm, n, :qp_pump, a), var(wm, n, :g_pump, a), var(wm, n, :z_pump, a)

    # Define the (relaxed) head gain relationship for the pump.
    c_1 = JuMP.@constraint(wm.model, g <= pc[1]*qp^2 + pc[2]*qp + pc[3]*z)

    # Add a linear lower bound on the head gain approximation.
    qp_ub = JuMP.upper_bound(qp)
    g_1, g_2 = pc[3], pc[1]*qp_ub^2 + pc[2]*qp_ub + pc[3]
    g_lb_line = (g_2 - g_1) * inv(qp_ub) * qp + g_1 * z
    c_2 = JuMP.@constraint(wm.model, g_lb_line <= g)

    append!(con(wm, n, :on_off_pump_head_gain, a), [c_1, c_2])
end


function objective_wf(wm::CQRDWaterModel)
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end


function objective_owf(wm::CQRDWaterModel)
    # Initialize the objective function.
    objective = JuMP.AffExpr(0.0)

    for (n, nw_ref) in nws(wm)
        # Get common constant parameters.
        rho = 1000.0 # Water density (kilogram per cubic meter).
        gravity = 9.80665 # Gravitational acceleration (meter per second squared).
        constant = rho * gravity * ref(wm, n, :time_step)

        for (a, pump) in nw_ref[:pump]
            q_min_active = get(pump, "q_min_active", _q_eps)

            if haskey(pump, "energy_price")
                # Get price and pump curve data.
                price = pump["energy_price"]
                head_curve = ref(wm, n, :pump, a)["head_curve"]
                curve_fun = _get_function_from_head_curve(head_curve)

                # Get flow-related variables and data.
                z = var(wm, n, :z_pump, a)
                qp, g = var(wm, n, :qp_pump, a), var(wm, n, :g_pump, a)
                points = collect(range(q_min_active, stop=JuMP.upper_bound(qp), length=50))

                # Get pump efficiency data.
                if haskey(pump, "efficiency_curve")
                    # Use the efficiency at the midpoint of the pump flow bounds.
                    eff_curve = pump["efficiency_curve"]
                    eff = _calc_efficiencies(points, eff_curve)
                else
                    eff = pump["efficiency"]
                end

                # Compute discrete costs from existing points.
                flows_cubed = _calc_cubic_flow_values(points, curve_fun)
                costs = (constant*price) .* inv.(eff) .* flows_cubed

                # Fit a quadratic function to the above discrete costs.
                LsqFit.@. func(x, p) = p[1]*x + p[2]
                fit = LsqFit.curve_fit(func, points, costs, zeros(length(costs)))
                coeffs = LsqFit.coef(fit)

                # Add the cost corresponding to the current pump's operation.
                JuMP.add_to_expression!(objective, coeffs[1]*qp + coeffs[2]*z)
            else
                Memento.error(_LOGGER, "No cost given for pump \"$(pump["name"])\"")
            end
        end
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end
