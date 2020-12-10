# Define common QRD (quadratic relaxation- and direction-based) implementations
# of water distribution network constraints, which use directed flow variables.


function _calc_pipe_exponent_second_order_expansion(q::Float64, exponent::Float64, v::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr})
    affine_term = q^exponent * z + exponent * q^(exponent - 1.0) * (v - q * z)
    quad_term = 0.5 * exponent * (exponent - 1.0) * q^(exponent - 2.0) * (v^2 + q^2 * z - 2.0 * v * q)
    return q != 0.0 ? affine_term + quad_term : 0.0
end


function constraint_pipe_head_loss(
    wm::AbstractQRDModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Get the flow direction variable.
    y = var(wm, n, :y_pipe, a)

    # Get variables for positive flow and head difference.
    qp, dhp = var(wm, n, :qp_pipe, a), var(wm, n, :dhp_pipe, a)

    # Get a second-order expansion of the positive nonlinearity at a point.
    qp_mid = 0.5 * JuMP.upper_bound(qp) # Get the positive flow point.
    expansion_p = _calc_pipe_exponent_second_order_expansion(qp_mid, exponent, qp, y)

    # Add constraints for the positive head loss approximation.
    c_1 = JuMP.@constraint(wm.model, r * expansion_p <= inv(L) * dhp)
    c_2 = JuMP.@constraint(wm.model, r * expansion_p >= inv(L) * dhp)

    # Get variables for negative flow and head difference.
    qn, dhn = var(wm, n, :qn_pipe, a), var(wm, n, :dhn_pipe, a)

    # Get a second-order expansion of the head loss nonlinearity at a point.
    qn_mid = 0.5 * JuMP.upper_bound(qn) # Get the negative flow point.
    expansion_n = _calc_pipe_exponent_second_order_expansion(qn_mid, exponent, qn, 1.0 - y)

    # Add constraints for the negative head loss approximation.
    c_3 = JuMP.@constraint(wm.model, r * expansion_n <= inv(L) * dhn)
    c_4 = JuMP.@constraint(wm.model, r * expansion_n >= inv(L) * dhn)

    # Append the :pipe_head_loss constraint array.
    append!(con(wm, n, :pipe_head_loss)[a], [c_1, c_2, c_3, c_4])
end


function objective_owf_default(wm::QRDWaterModel)
    # Initialize the objective function.
    objective = zero(JuMP.QuadExpr)

    for (n, nw_ref) in nws(wm)
        # Get common constant parameters.
        constant = _DENSITY * _GRAVITY * ref(wm, n, :time_step)

        for (a, pump) in nw_ref[:pump]
            q_min_forward = get(pump, "flow_min_forward", _FLOW_MIN)

            if haskey(pump, "energy_price")
                # Get price and pump curve data.
                price = pump["energy_price"]
                head_curve = ref(wm, n, :pump, a)["head_curve"]
                curve_fun = _get_function_from_head_curve(head_curve)

                # Get flow-related variables and data.
                z = var(wm, n, :z_pump, a)
                qp, g = var(wm, n, :qp_pump, a), var(wm, n, :g_pump, a)
                points = collect(range(q_min_forward, stop=JuMP.upper_bound(qp), length=50))

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
                LsqFit.@. func(x, p) = p[1]*x*x + p[2]*x + p[3]
                fit = LsqFit.curve_fit(func, points, costs, zeros(length(costs)))
                coeffs = LsqFit.coef(fit)

                # Add the cost corresponding to the current pump's operation.
                JuMP.add_to_expression!(objective, coeffs[1], qp, qp)
                JuMP.add_to_expression!(objective, coeffs[2] * qp + coeffs[3] * z)
            else
                Memento.error(_LOGGER, "No cost given for pump \"$(pump["name"])\"")
            end
        end
    end

    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end
