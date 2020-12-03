function _calc_pump_flow_min(pump::Dict{String,<:Any}, node_fr::Dict{String,Any}, node_to::Dict{String,Any})
    return max(0.0, get(pump, "flow_min", 0.0))
end


function _calc_pump_flow_min_forward(pump::Dict{String,<:Any}, node_fr::Dict{String,Any}, node_to::Dict{String,Any})
    flow_min_forward = get(pump, "flow_min_forward", _FLOW_MIN)
    return max(_calc_pump_flow_min(pump, node_fr, node_to), flow_min_forward)
end


function _calc_pump_flow_max_reverse(pump::Dict{String,<:Any}, node_fr::Dict{String,Any}, node_to::Dict{String,Any})
    flow_max_reverse = get(pump, "flow_max_reverse", 0.0)
    return min(_calc_pump_flow_max(pump, node_fr, node_to), flow_max_reverse)
end


function _calc_pump_flow_max(pump::Dict{String,<:Any}, node_fr::Dict{String,Any}, node_to::Dict{String,Any})
    coeff = _get_function_from_head_curve(pump["head_curve"])
    q_max_1 = (-coeff[2] + sqrt(coeff[2]^2 - 4.0*coeff[1]*coeff[3])) * inv(2.0*coeff[1])
    q_max_2 = (-coeff[2] - sqrt(coeff[2]^2 - 4.0*coeff[1]*coeff[3])) * inv(2.0*coeff[1])
    return min(max(q_max_1, q_max_2), get(pump, "flow_max", Inf))
end


function _calc_pump_flow_bounds_active(pump::Dict{String,<:Any})
    q_min, q_max = _calc_pump_flow_bounds(pump)
    q_min = max(max(get(pump, "flow_min_forward", _FLOW_MIN), _FLOW_MIN), q_min)
    return q_min, q_max
end


function _calc_pump_best_efficiency_curve(pump::Dict{String, <:Any})
    # Build a two-dimensional array of the efficiency curve points.
    efficiency_array = vcat([hcat(x[1], x[2]) for x in pump["efficiency_curve"]]...)

    # Build another array for fitting the efficiency curve.
    fit_array = hcat(efficiency_array[:, 1].^2, efficiency_array[:, 1])

    # Perform a fit of the efficiency curve and get the linear coefficients.
    return fit_array \ efficiency_array[:, 2]
end


function _calc_pump_best_efficiency_head_curve(pump::Dict{String, <:Any})
    # Build a two-dimensional array of the head curve points.
    head_array = vcat([hcat(x[1], x[2]) for x in pump["head_curve"]]...)

    # Build another array for fitting the efficiency curve.
    fit_array = hcat(head_array[:, 1].^2, ones(size(head_array, 1)))

    # Perform a fit of the head curve and get the linear coefficients.
    return fit_array \ head_array[:, 2]
end


function _calc_pump_best_efficiency(pump::Dict{String, <:Any})
    if haskey(pump, "efficiency_curve")
        # Perform a fit of the efficiency curve and get the linear coefficients.
        coeffs = _calc_pump_best_efficiency_curve(pump)

        # Return the maximum efficiency predicted by the fitted curve.
        return -0.25 * coeffs[2]^2 * inv(coeffs[1])
    else
        # Assume the single value specified for efficiency is the best efficiency.
        return pump["efficiency"]
    end
end


function _calc_pump_best_efficiency_flow(pump::Dict{String, <:Any})
    if haskey(pump, "efficiency_curve")
        # Perform a fit of the efficiency curve and get the linear coefficients.
        coeffs = _calc_pump_best_efficiency_curve(pump)

        # Return the flow corresponding to the best efficiency point on the curve.
        return -0.5*coeffs[2] * inv(coeffs[1])
    else
        # An efficiency curve was not provided. Flow must be determined from the head curve.
        if length(pump["head_curve"]) == 1
            return pump["head_curve"][1][1]
        else
            # Perform a fit of the head curve and get the linear coefficients.
            coeffs = _calc_pump_best_efficiency_head_curve(pump)

            # Return the flow at which the maximum head gain occurs.
            return sqrt(-0.25 * coeffs[2] * inv(coeffs[1]))
        end
    end
end


function _calc_pump_best_efficiency_head_gain(pump::Dict{String, <:Any})
    if haskey(pump, "efficiency_curve")
        # We could determine the best efficiency flow, q, from the head curve and check
        # whether it agrees with the value determined from the efficiency curve. However, we
        # will assume that the one from the efficiency curve is more accurate.
        q = _calc_pump_best_efficiency_flow(pump)

        # Build a two-dimensional array of the head curve points.
        head_array = vcat([hcat(x[1], x[2]) for x in pump["head_curve"]]...)

        # Build another array for fitting the efficiency curve.
        fit_array = -inv(3.0) * inv(q^2) * head_array[:, 1].^2 .+ (4.0 * inv(3.0))

        # Return the head gain predicted by the best efficiency curve.
        return fit_array \ head_array[:, 2]
    else
        # An efficiency curve was not provided. Gain must be determined from the head curve.
        if length(pump["head_curve"]) == 1
            return pump["head_curve"][1][2]
        else
            # Perform a fit of the head curve and get the linear coefficients.
            coeffs = _calc_pump_best_efficiency_head_curve(pump)

            # Return the head at which the maximum head gain occurs.
            return 0.75 * coeffs[2]
        end
    end
end


function _calc_pump_best_efficiency_power(pump::Dict{String, <:Any})
    efficiency = _calc_pump_best_efficiency(pump)
    flow = _calc_pump_best_efficiency_flow(pump)
    head_gain = _calc_pump_best_efficiency_head_gain(pump)
    return _DENSITY * _GRAVITY * inv(efficiency) * flow * head_gain
end


function _calc_pump_best_efficiency_head_gain_curve(pump::Dict{String, <:Any})
    flow = _calc_pump_best_efficiency_flow(pump)
    head_gain = _calc_pump_best_efficiency_head_gain(pump)
    return [-inv(3.0) * head_gain * inv(flow^2), 0.0, 4.0 * head_gain * inv(3.0)]
end


function _calc_pump_head_gain_curve(pump::Dict{String, <:Any})
    head_curve = pump["head_curve"]

    LsqFit.@. func(x, p) = p[1]*x*x + p[2]*x + p[3]

    if length(head_curve) > 1
        fit = LsqFit.curve_fit(func, first.(head_curve), last.(head_curve), [0.0, 0.0, 0.0])
        return LsqFit.coef(fit)
    elseif length(head_curve) == 1
        new_points = [(0.0, 1.33 * head_curve[1][2]), (2.0 * head_curve[1][1], 0.0)]
        head_curve = vcat(new_points, head_curve)
        fit = LsqFit.curve_fit(func, first.(head_curve), last.(head_curve), [0.0, 0.0, 0.0])
        return LsqFit.coef(fit)
    end
end
