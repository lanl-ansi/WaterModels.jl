function correct_pumps!(data::Dict{String, <:Any})
    for (idx, pump) in data["pump"]
        # Get common connecting node data for later use.
        node_fr = data["node"][string(pump["node_fr"])]
        node_to = data["node"][string(pump["node_to"])]

        # Correct various pump properties. The sequence is important, here.
        _correct_pump_head_curve_form!(pump)
        _correct_pump_flow_bounds!(pump, node_fr, node_to)
    end
end


function _correct_pump_head_curve_form!(pump::Dict{String, <:Any})
    pump["head_curve_form"] = get(pump, "head_curve_form", QUADRATIC)
end


function _correct_pump_flow_bounds!(pump::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any})
    pump["flow_min"] = _calc_pump_flow_min(pump, node_fr, node_to)
    pump["flow_max"] = _calc_pump_flow_max(pump, node_fr, node_to)
    pump["flow_min_forward"] = _calc_pump_flow_min_forward(pump, node_fr, node_to)
    pump["flow_max_reverse"] = _calc_pump_flow_max_reverse(pump, node_fr, node_to)
end


function _calc_pump_flow_min(pump::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any})
    return max(0.0, get(pump, "flow_min", 0.0))
end


function _calc_pump_flow_min_forward(pump::Dict{String, <:Any}, node_fr::Dict{String,<: Any}, node_to::Dict{String, <:Any})
    flow_min_forward = get(pump, "flow_min_forward", _FLOW_MIN)
    return max(_calc_pump_flow_min(pump, node_fr, node_to), flow_min_forward)
end


function _calc_pump_flow_max_reverse(pump::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any})
    flow_max_reverse = get(pump, "flow_max_reverse", 0.0)
    return min(_calc_pump_flow_max(pump, node_fr, node_to), flow_max_reverse)
end


function calc_pump_head_gain_max(pump::Dict{String,<:Any}, node_fr::Dict{String,Any}, node_to::Dict{String,Any})
    # Calculate the flow at the maximum head gain, then return maximum head gain.
    c = _calc_head_curve_coefficients(pump)
    flow_at_max = -c[2] * inv(2.0 * c[1]) > 0.0 ? -c[2] * inv(2.0 * c[1]) : 0.0
    return c[1]*flow_at_max^2 + c[2]*flow_at_max + c[3]
end


function _calc_pump_flow_max(pump::Dict{String,<:Any}, node_fr::Dict{String,Any}, node_to::Dict{String,Any})
    # Get possible maximal flow values based on the head curve.
    c = _calc_head_curve_coefficients(pump)
    q_max_1 = (-c[2] + sqrt(c[2]^2 - 4.0*c[1]*c[3])) * inv(2.0*c[1])
    q_max_2 = (-c[2] - sqrt(c[2]^2 - 4.0*c[1]*c[3])) * inv(2.0*c[1])

    # Get possible maximal flow values based on maximum head gain.
    g = get(node_to, "head_max", Inf) - get(node_fr, "head_min", -Inf)
    q_max_3 = g < Inf ? (-c[2] - sqrt(c[2]^2 - 4.0*c[1]*(c[3] + g))) * inv(2.0*c[1]) : Inf
    q_max_4 = g < Inf ? (-c[2] + sqrt(c[2]^2 - 4.0*c[1]*(c[3] + g))) * inv(2.0*c[1]) : Inf

    # Get the minimal value of the above and the possible "flow_max" value.
    return min(max(q_max_1, q_max_2), max(q_max_3, q_max_4), get(pump, "flow_max", Inf))
end


function _calc_head_curve_coefficients(pump::Dict{String, <:Any})
    if pump["head_curve_form"] == QUADRATIC
        return _calc_head_curve_coefficients_quadratic(pump)
    elseif pump["head_curve_form"] == BEST_EFFICIENCY_POINT
        return _calc_head_curve_coefficients_best_efficiency_point(pump)
    else
        error("\"$(pump["head_curve_form"])\" is not a valid head curve formulation.")
    end
end


function _calc_head_curve_function(pump::Dict{String, <:Any})
    if pump["head_curve_form"] == QUADRATIC
        coeff = _calc_head_curve_coefficients_quadratic(pump)
        return x -> coeff .* [x^2, x, 1.0]
    elseif pump["head_curve_form"] == BEST_EFFICIENCY_POINT
        coeff = _calc_head_curve_coefficients_best_efficiency_point(pump)
        return x -> coeff .* [x^2, x, 1.0]
    else
        error("\"$(pump["head_curve_form"])\" is not a valid head curve formulation.")
    end
end


function _calc_head_curve_coefficients_quadratic(pump::Dict{String, <:Any})
    if length(pump["head_curve"]) > 1
        array = pump["head_curve"]
    elseif length(pump["head_curve"]) == 1
        array = [0.0 1.33 * pump["head_curve"][1][2]; 2.0 * pump["head_curve"][1][1] 0.0]
    else
        error("Pump \"$(pump["name"])\" has no head curve points.")
    end

    # Build a two-dimensional array of the head curve points.
    array = vcat([hcat(x[1], x[2]) for x in pump["head_curve"]]...)

    # Build another array for fitting the efficiency curve.
    fit_array = hcat(array[:, 1].^2, array[:, 1], ones(size(array, 1)))

    # Perform a fit of the head curve and return the model coefficients.
    return fit_array \ array[:, 2]
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


function _calc_head_curve_coefficients_best_efficiency_point(pump::Dict{String, <:Any})
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
