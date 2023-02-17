function correct_pumps!(data::Dict{String, <:Any})
    wm_data = get_wm_data(data)

    if wm_data["per_unit"]
        flow_transform = _calc_flow_per_unit_transform(wm_data)
        flow_epsilon = flow_transform(_FLOW_MIN)
        func = x -> _correct_pumps!(x, flow_epsilon)
    else
        func = x -> _correct_pumps!(x, _FLOW_MIN)
    end

    apply_wm!(func, data; apply_to_subnetworks = true)
end

function correct_ne_pumps!(data::Dict{String, <:Any})
    wm_data = get_wm_data(data)

    if wm_data["per_unit"]
        flow_transform = _calc_flow_per_unit_transform(wm_data)
        flow_epsilon = flow_transform(_FLOW_MIN)
        func = x -> _correct_ne_pumps!(x, flow_epsilon)
    else
        func = x -> _correct_ne_pumps!(x, _FLOW_MIN)
    end

    apply_wm!(func, data; apply_to_subnetworks = true)
end

function _correct_pumps!(data::Dict{String, <:Any}, flow_epsilon::Float64)
    for (idx, pump) in data["pump"]
        # Get common connecting node data for later use.
        node_fr = data["node"][string(pump["node_fr"])]
        node_to = data["node"][string(pump["node_to"])]

        # Correct various pump properties. The sequence is important, here.
        _correct_status!(pump)
        _correct_flow_direction!(pump)
        _correct_pump_type!(pump)
        _correct_pump_flow_bounds!(pump, node_fr, node_to, flow_epsilon)
    end
end

function _correct_ne_pumps!(data::Dict{String, <:Any}, flow_epsilon::Float64)
    for (idx, pump) in data["ne_pump"]
        # Get common connecting node data for later use.
        node_fr = data["node"][string(pump["node_fr"])]
        node_to = data["node"][string(pump["node_to"])]

        # Correct various expansion pump properties. The sequence is important, here.
        _correct_ne_pump_data!(pump)
        _correct_status!(pump)
        _correct_flow_direction!(pump)
        _correct_pump_type!(pump)
        _correct_pump_flow_bounds!(pump, node_fr, node_to, flow_epsilon)
    end
end

function _correct_ne_pump_data!(ne_pump:: Dict{String, <:Any})
    ne_pump["head_curve"] = eval(Meta.parse(ne_pump["head_curve"]))
    ne_pump["pump_type"] = eval(Meta.parse(ne_pump["pump_type"]))
    ne_pump["status"] = eval(Meta.parse(ne_pump["status"]))
    ne_pump["flow_direction"] = eval(Meta.parse(ne_pump["flow_direction"]))
end

function set_pump_flow_partition!(
    pump::Dict{String, <:Any}, error_tolerance::Float64, length_tolerance::Float64)
    # Compute the head gain function and its derivative.
    f = _calc_head_curve_function(pump)
    f_dash = _calc_head_curve_derivative(pump)

    # Initialize the partitioning of flows for the pipe.
    partition = Vector{Float64}([pump["flow_min_forward"], pump["flow_max"]])

    # Use PolyhedralRelaxations to determine partitions with desired accuracy.
    uvf_data = PolyhedralRelaxations.UnivariateFunctionData(
        f, f_dash, partition, error_tolerance,
        length_tolerance, 1.0e-6, 9e9, length(partition))
    PolyhedralRelaxations._refine_partition!(uvf_data)

    # Set pump flow partition using the above partitioning.
    pump["flow_partition"] = Vector{Float64}(uvf_data.partition)
end


function correct_pump_types!(data::Dict{String,<:Any})
    apply_wm!(_correct_pump_types!, data; apply_to_subnetworks = true)
end


function correct_ne_pump_types!(data::Dict{String,<:Any})
    apply_wm!(_correct_ne_pump_types!, data; apply_to_subnetworks = true)
end


function _correct_pump_types!(data::Dict{String,<:Any})
    components = values(get(data, "pump", Dict{String,Any}()))
    _correct_pump_type!.(components)
end

function _correct_ne_pump_types!(data::Dict{String,<:Any})
    components = values(get(data, "ne_pump", Dict{String,Any}()))
    _correct_ne_pump_type!.(components)
end


function _correct_pump_type!(pump::Dict{String, <:Any})
    head_curve_tmp = get(pump, "pump_type", PUMP_QUADRATIC)

    if isa(head_curve_tmp, PUMP)
        pump["pump_type"] = head_curve_tmp
    else
        pump["pump_type"] = PUMP(head_curve_tmp)
    end
end

function _correct_ne_pump_type!(ne_pump::Dict{String, <:Any})
    head_curve_tmp = get(ne_pump, "pump_type", PUMP_QUADRATIC)

    if isa(head_curve_tmp, PUMP)
        ne_pump["pump_type"] = head_curve_tmp
    else
        ne_pump["pump_type"] = PUMP(head_curve_tmp)
    end
end


function _correct_pump_flow_bounds!(pump::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any}, flow_epsilon::Float64)
    if get(pump, "z_max", 1.0) == 0.0
        pump["flow_min"] = 0.0
        pump["flow_max"] = 0.0
        pump["flow_min_forward"] = 0.0
        pump["flow_max_reverse"] = 0.0
    else
        pump["flow_min"] = max(0.0, _calc_pump_flow_min(pump, node_fr, node_to))
        pump["flow_max"] = max(0.0, _calc_pump_flow_max(pump, node_fr, node_to))
        pump["flow_min_forward"] = max(0.0, max(pump["flow_min"], _calc_pump_flow_min_forward(pump, node_fr, node_to, flow_epsilon)))
        pump["flow_min_forward"] = min(pump["flow_max"], pump["flow_min_forward"])
        pump["flow_max_reverse"] = max(0.0, min(pump["flow_max"], _calc_pump_flow_max_reverse(pump, node_fr, node_to)))
    end

    pump["flow_max"] = max(pump["flow_min"], pump["flow_max"])
    pump["flow_min"] = min(pump["flow_min"], pump["flow_max"])
    pump["flow_min_forward"] = max(pump["flow_min_forward"], pump["flow_min"])
    pump["flow_min_forward"] = min(pump["flow_min_forward"], pump["flow_max"])
    pump["flow_max_reverse"] = min(pump["flow_max_reverse"], pump["flow_max"])
    pump["flow_max_reverse"] = max(pump["flow_max_reverse"], pump["flow_min"])

    @assert pump["flow_min"] <= pump["flow_max"]
    @assert get(pump, "flow_min_forward", 0.0) <= max(0.0, pump["flow_max"])
    @assert min(0.0, pump["flow_min"]) <= get(pump, "flow_max_reverse", 0.0)
end


function _calc_pump_flow_min(pump::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any})
    return max(0.0, get(pump, "flow_min", 0.0))
end


function _calc_pump_flow_min_forward(pump::Dict{String, <:Any}, node_fr::Dict{String,<: Any}, node_to::Dict{String, <:Any}, flow_epsilon::Float64)
    flow_min_forward = get(pump, "flow_min_forward", flow_epsilon)
    return max(_calc_pump_flow_min(pump, node_fr, node_to), flow_min_forward)
end


function _calc_pump_flow_max_reverse(pump::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any})
    flow_max_reverse = get(pump, "flow_max_reverse", 0.0)
    return min(_calc_pump_flow_max(pump, node_fr, node_to), flow_max_reverse)
end


function _calc_pump_head_gain_max(pump::Dict{String, <:Any}, node_fr::Dict{String,Any}, node_to::Dict{String,Any})
    # Calculate the flow at the maximum head gain, then return maximum head gain.
    c = _calc_head_curve_coefficients(pump)

    if pump["pump_type"] in [PUMP_QUADRATIC, PUMP_BEST_EFFICIENCY_POINT]
        flow_at_max = -c[2] * inv(2.0 * c[1]) > 0.0 ? -c[2] * inv(2.0 * c[1]) : 0.0
        return max(0.0, c[1] * flow_at_max^2 + c[2] * flow_at_max + c[3])
    elseif pump["pump_type"] in [PUMP_EPANET, PUMP_LINEAR_POWER]
        return max(0.0, c[1])
    end
end


function _calc_pump_flow_max(pump::Dict{String,<:Any}, node_fr::Dict{String,Any}, node_to::Dict{String,Any})
    # Get possible maximal flow values based on the head curve.
    c = _calc_head_curve_coefficients(pump)

    if pump["pump_type"] in [PUMP_QUADRATIC, PUMP_BEST_EFFICIENCY_POINT]
        q_max_1 = (-c[2] + sqrt(c[2]^2 - 4.0*c[1]*c[3])) * inv(2.0 * c[1])
        q_max_2 = (-c[2] - sqrt(c[2]^2 - 4.0*c[1]*c[3])) * inv(2.0 * c[1])

        # Get possible maximal flow values based on maximum head gain.
        g = get(node_to, "head_max", Inf) - get(node_fr, "head_min", -Inf)
        q_max_3 = g < Inf ? (-c[2] - sqrt(c[2]^2 - 4.0*c[1]*(c[3] + g))) * inv(2.0*c[1]) : Inf
        q_max_4 = g < Inf ? (-c[2] + sqrt(c[2]^2 - 4.0*c[1]*(c[3] + g))) * inv(2.0*c[1]) : Inf

        # Get the minimal value of the above and the possible "flow_max" value.
        return min(max(q_max_1, q_max_2), max(q_max_3, q_max_4), get(pump, "flow_max", Inf))
    elseif pump["pump_type"] in [PUMP_EPANET, PUMP_LINEAR_POWER]
        return min((-c[1] * inv(c[2]))^(inv(c[3])), get(pump, "flow_max", Inf))
    end
end


function _calc_pump_power_max(
    pump::Dict{String,<:Any}, node_fr::Dict{String,Any}, node_to::Dict{String,Any},
    density_scaled::Float64, gravity_scaled::Float64)
    flow_max = _calc_pump_flow_max(pump, node_fr, node_to)
    gain_max = _calc_pump_head_gain_max(pump, node_fr, node_to)

    if haskey(pump, "efficiency_curve")
        min_efficiency = minimum(x[2] for x in pump["efficiency_curve"])
    else
        min_efficiency = pump["efficiency"]
    end

    # Ensure minimum efficiency is greater than zero.
    @assert min_efficiency > 0.0

    # Return the maximum pump power.
    return density_scaled * gravity_scaled *
        flow_max * gain_max * inv(min_efficiency)
end


function _calc_head_curve_coefficients(pump::Dict{String, <:Any})
    if pump["pump_type"] == PUMP_QUADRATIC
        return _calc_head_curve_coefficients_quadratic(pump)
    elseif pump["pump_type"] == PUMP_BEST_EFFICIENCY_POINT
        return _calc_head_curve_coefficients_best_efficiency_point(pump)
    elseif pump["pump_type"] in [PUMP_EPANET, PUMP_LINEAR_POWER]
        return _calc_head_curve_coefficients_epanet(pump)
    else
        error("\"$(pump["pump_type"])\" is not a valid head curve formulation.")
    end
end


function _calc_head_curve_function(pump::Dict{String, <:Any})
    if pump["pump_type"] == PUMP_QUADRATIC
        coeff = _calc_head_curve_coefficients_quadratic(pump)
        return x -> sum(coeff .* [x^2, x, 1.0])
    elseif pump["pump_type"] == PUMP_BEST_EFFICIENCY_POINT
        coeff = _calc_head_curve_coefficients_best_efficiency_point(pump)
        return x -> sum(coeff .* [x^2, x, 1.0])
    elseif pump["pump_type"] in [PUMP_EPANET, PUMP_LINEAR_POWER]
        coeff = _calc_head_curve_coefficients_epanet(pump)
        return x -> coeff[1] + coeff[2] * x^coeff[3]
    else
        error("\"$(pump["pump_type"])\" is not a valid head curve formulation.")
    end
end


function _calc_head_curve_function(pump::Dict{String, <:Any}, z::JuMP.VariableRef)
    if pump["pump_type"] == PUMP_QUADRATIC
        coeff = _calc_head_curve_coefficients_quadratic(pump)
        return x -> sum(coeff .* [x^2, x, z])
    elseif pump["pump_type"] == PUMP_BEST_EFFICIENCY_POINT
        coeff = _calc_head_curve_coefficients_best_efficiency_point(pump)
        return x -> sum(coeff .* [x^2, x, z])
    elseif pump["pump_type"] in [PUMP_EPANET, PUMP_LINEAR_POWER]
        coeff = _calc_head_curve_coefficients_epanet(pump)
        return x -> coeff[1] * z + coeff[2] * x^coeff[3]
    else
        error("\"$(pump["pump_type"])\" is not a valid head curve formulation.")
    end
end

function _calc_head_curve_derivative(pump::Dict{String, <:Any})
    if pump["pump_type"] == PUMP_QUADRATIC
        coeff = _calc_head_curve_coefficients_quadratic(pump)
        return x -> sum(coeff .* [2.0 * x, 1.0, 0.0])
    elseif pump["pump_type"] == PUMP_BEST_EFFICIENCY_POINT
        coeff = _calc_head_curve_coefficients_best_efficiency_point(pump)
        return x -> sum(coeff .* [2.0 * x, 1.0, 0.0])
    elseif pump["pump_type"] in [PUMP_EPANET, PUMP_LINEAR_POWER]
        coeff = _calc_head_curve_coefficients_epanet(pump)
        return x -> coeff[2] * coeff[3] * x^(coeff[3] - 1.0)
    else
        error("\"$(pump["pump_type"])\" is not a valid head curve formulation.")
    end
end


function _calc_head_curve_coefficients_epanet(pump::Dict{String, <:Any})
    # Use LsqFit to compute a fit with respect to points on the head curve.
    q, h = [x[1] for x in pump["head_curve"]], [x[2] for x in pump["head_curve"]]
    model_function(x, p) = p[1] .+ p[2] .* x.^p[3] # (i.e., a + b * x^c)
    params = LsqFit.curve_fit(model_function, q, h, [h[1], -1.0, 2.0];
        lower = [0.0, -Inf, 1.0 + 1.0e-12], upper = [Inf, -1.0e-12, Inf]).param

    # Ensure the function is concave with a positive offset.
    @assert params[1] > 0.0 && params[2] * (params[3]^2 - params[3]) < 0.0

    # Return the vector of function parameters.
    return Vector{Float64}(params)
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

    # Build another array for fitting the head curve.
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
        return -0.5 * coeffs[2] * inv(coeffs[1])
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


function _calc_pump_best_efficiency_power(pump::Dict{String, <:Any}, density::Float64, gravity::Float64)
    efficiency = _calc_pump_best_efficiency(pump)
    flow = _calc_pump_best_efficiency_flow(pump)
    head_gain = _calc_pump_best_efficiency_head_gain(pump)
    return density * gravity * inv(efficiency) * flow * head_gain
end


function _calc_head_curve_coefficients_best_efficiency_point(pump::Dict{String, <:Any})
    flow = _calc_pump_best_efficiency_flow(pump)
    head_gain = _calc_pump_best_efficiency_head_gain(pump)
    return [-inv(3.0) * head_gain * inv(flow^2), 0.0, 4.0 * head_gain * inv(3.0)]
end


function _calc_pump_power_points(wm::AbstractWaterModel, nw::Int, pump_id::Int, num_points::Int)
    pump = ref(wm, nw, :pump, pump_id)
    head_curve_function = ref(wm, nw, :pump, pump_id, "head_curve_function")

    wm_data = get_wm_data(wm.data)
    flow_transform = _calc_flow_per_unit_transform(wm_data)
    q_build = range(0.0, stop = pump["flow_max"] + 1.0e-7, length = num_points)
    f_build = head_curve_function.(collect(q_build)) .* q_build

    if haskey(pump, "efficiency_curve")
        eff_curve = pump["efficiency_curve"]
        eff = _calc_efficiencies(collect(q_build), eff_curve)
    else
        eff = pump["efficiency"]
    end

    base_mass = 1.0 / _calc_mass_per_unit_transform(wm_data)(1.0)
    base_time = 1.0 / _calc_time_per_unit_transform(wm_data)(1.0)
    base_length = 1.0 / _calc_length_per_unit_transform(wm_data)(1.0)
    flow_transform = 1.0 / _calc_flow_per_unit_transform(wm_data)(1.0)

    density = _calc_scaled_density(base_mass, base_length)
    gravity = _calc_scaled_gravity(base_length, base_time)
    return q_build, max.(0.0, density * gravity * inv.(eff) .* f_build)
end

function _calc_pump_power_points_ne(wm::AbstractWaterModel, nw::Int, pump_id::Int, num_points::Int)
    pump = ref(wm, nw, :ne_pump, pump_id)
    head_curve_function = ref(wm, nw, :ne_pump, pump_id, "head_curve_function")

    wm_data = get_wm_data(wm.data)
    flow_transform = _calc_flow_per_unit_transform(wm_data)
    q_build = range(0.0, stop = pump["flow_max"] + 1.0e-7, length = num_points)
    f_build = head_curve_function.(collect(q_build)) .* q_build

    if haskey(pump, "efficiency_curve")
        eff_curve = pump["efficiency_curve"]
        eff = _calc_efficiencies(collect(q_build), eff_curve)
    else
        eff = pump["efficiency"]
    end

    base_mass = 1.0 / _calc_mass_per_unit_transform(wm_data)(1.0)
    base_time = 1.0 / _calc_time_per_unit_transform(wm_data)(1.0)
    base_length = 1.0 / _calc_length_per_unit_transform(wm_data)(1.0)
    flow_transform = 1.0 / _calc_flow_per_unit_transform(wm_data)(1.0)

    density = _calc_scaled_density(base_mass, base_length)
    gravity = _calc_scaled_gravity(base_length, base_time)
    return q_build, max.(0.0, density * gravity * inv.(eff) .* f_build)
end


function _calc_pump_power(wm::AbstractWaterModel, nw::Int, pump_id::Int, q::Vector{Float64})
    q_true, f_true = _calc_pump_power_points(wm, nw, pump_id, 100)
    return max.(Interpolations.LinearInterpolation(q_true, f_true).(q), 0.0)
end

function _calc_pump_power_ne(wm::AbstractWaterModel, nw::Int, pump_id::Int, q::Vector{Float64})
    q_true, f_true = _calc_pump_power_points_ne(wm, nw, pump_id, 100)
    return max.(Interpolations.LinearInterpolation(q_true, f_true).(q), 0.0)
end

function _calc_pump_power_ua(wm::AbstractWaterModel, nw::Int, pump_id::Int, q::Vector{Float64})
    q_true, f_true = _calc_pump_power_points(wm, nw, pump_id, 100)
    f_interp = Interpolations.LinearInterpolation(q_true, f_true).(q)

    for i in 2:length(q)
        # Find the indices that will be used in the approximation.
        true_ids = filter(k -> q_true[k] >= q[i-1] && q_true[k] <= q[i], 1:length(q_true))

        if length(true_ids) == 0
            true_ids = filter(k -> q_true[k] >= q[i-1], 1:length(q_true))
            true_ids = vcat(true_ids[1] - 1, true_ids[1])
        end

        # Shift the approximation by the maximum error, ensuring it's an underapproximation.
        slope = (f_interp[i] - f_interp[i-1]) / (q[i] - q[i-1])
        f_est_s = f_interp[i-1] .+ (slope .* (q_true[true_ids] .- q[i-1]))
        est_err = max(0.0, maximum(f_est_s .- f_true[true_ids]))
        f_interp[i-1:i] .-= est_err
    end

    return f_interp
end

function _calc_pump_power_ua_ne(wm::AbstractWaterModel, nw::Int, pump_id::Int, q::Vector{Float64})
    q_true, f_true = _calc_pump_power_points_ne(wm, nw, pump_id, 100)
    f_interp = Interpolations.LinearInterpolation(q_true, f_true).(q)

    for i in 2:length(q)
        # Find the indices that will be used in the approximation.
        true_ids = filter(k -> q_true[k] >= q[i-1] && q_true[k] <= q[i], 1:length(q_true))

        if length(true_ids) == 0
            true_ids = filter(k -> q_true[k] >= q[i-1], 1:length(q_true))
            true_ids = vcat(true_ids[1] - 1, true_ids[1])
        end

        # Shift the approximation by the maximum error, ensuring it's an underapproximation.
        slope = (f_interp[i] - f_interp[i-1]) / (q[i] - q[i-1])
        f_est_s = f_interp[i-1] .+ (slope .* (q_true[true_ids] .- q[i-1]))
        est_err = max(0.0, maximum(f_est_s .- f_true[true_ids]))
        f_interp[i-1:i] .-= est_err
    end

    return f_interp
end


function _calc_pump_power_oa(wm::AbstractWaterModel, nw::Int, pump_id::Int, q::Vector{Float64})
    q_true, f_true = _calc_pump_power_points(wm, nw, pump_id, 100)
    f_interp = Interpolations.LinearInterpolation(q_true, f_true).(q)

    for i in 2:length(q)
        slope = (f_interp[i] - f_interp[i-1]) / (q[i] - q[i-1])
        true_ids = filter(x -> q_true[x] >= q[i-1] && q_true[x] <= q[i], 1:length(q_true))
        f_est_s = f_interp[i-1] .+ (slope .* (q_true[true_ids] .- q[i-1]))
        est_err = max(0.0, maximum(f_true[true_ids] .- f_est_s))
        f_interp[i-1:i] .+= est_err
    end

    return f_interp
end

function _calc_pump_power_oa_ne(wm::AbstractWaterModel, nw::Int, pump_id::Int, q::Vector{Float64})
    q_true, f_true = _calc_pump_power_points_ne(wm, nw, pump_id, 100)
    f_interp = Interpolations.LinearInterpolation(q_true, f_true).(q)

    for i in 2:length(q)
        slope = (f_interp[i] - f_interp[i-1]) / (q[i] - q[i-1])
        true_ids = filter(x -> q_true[x] >= q[i-1] && q_true[x] <= q[i], 1:length(q_true))
        f_est_s = f_interp[i-1] .+ (slope .* (q_true[true_ids] .- q[i-1]))
        est_err = max(0.0, maximum(f_true[true_ids] .- f_est_s))
        f_interp[i-1:i] .+= est_err
    end

    return f_interp
end


function _calc_pump_power_quadratic_approximation(wm::AbstractWaterModel, nw::Int, pump_id::Int, z::JuMP.VariableRef)
    # Get good approximations of pump flow and power points.
    q_true, f_true = _calc_pump_power_points(wm, nw, pump_id, 100)

    # Build a two-dimensional array of the feature points.
    q_array = hcat(q_true .* q_true, q_true, ones(length(q_true)))

    # Obtain coefficients for the quadratic approximation from a least squares fit.
    linear_coefficients = q_array \ f_true

    # Ensure that a negative constant term does not exist.
    linear_coefficients[3] = max(0.0, linear_coefficients[3])

    # Return the least squares-fitted quadratic approximation.
    return x -> sum(linear_coefficients .* [x * x, x, z])
end

function _calc_pump_power_quadratic_approximation_ne(wm::AbstractWaterModel, nw::Int, pump_id::Int, z::JuMP.VariableRef)
    # Get good approximations of pump flow and power points.
    q_true, f_true = _calc_pump_power_points_ne(wm, nw, pump_id, 100)

    # Build a two-dimensional array of the feature points.
    q_array = hcat(q_true .* q_true, q_true, ones(length(q_true)))

    # Obtain coefficients for the quadratic approximation from a least squares fit.
    linear_coefficients = q_array \ f_true

    # Ensure that a negative constant term does not exist.
    linear_coefficients[3] = max(0.0, linear_coefficients[3])

    # Return the least squares-fitted quadratic approximation.
    return x -> sum(linear_coefficients .* [x * x, x, z])
end

function _calc_efficiencies(points::Vector{Float64}, curve::Vector{<:Any})
    q, eff = [[x[1] for x in curve], [x[2] for x in curve]]
    return Interpolations.LinearInterpolation(q, eff,
        extrapolation_bc=Interpolations.Flat()).(points)
end


function get_pump_flow_partition(pump::Dict{String, <:Any})
    @assert haskey(pump, "flow_partition")
    flows = filter(x -> x > 0.0, pump["flow_partition"])
    lower_bound = max(0.0, get(pump, "flow_min_forward", 0.0))

    if length(flows) > 0 && lower_bound == minimum(flows)
        return flows
    else
        flow_max = length(flows) > 0 ? maximum(flows) : lower_bound
        return lower_bound != flow_max ? vcat(lower_bound, flows) : [lower_bound]
    end
end


function get_pump_head_gain_partition(pump::Dict{String, <:Any})
    flow_partition = get_pump_flow_partition(pump)
    head_curve_function = _calc_head_curve_function(pump)
    return head_curve_function.(flow_partition)
end


function set_pump_warm_start!(data::Dict{String, <:Any})
    apply_wm!(_set_pump_warm_start!, data)
end

function set_ne_pump_warm_start!(data::Dict{String, <:Any})
    apply_wm!(_set_ne_pump_warm_start!, data)
end


function _relax_pumps!(data::Dict{String,<:Any})
    if !_IM.ismultinetwork(data)
        if haskey(data, "time_series") && haskey(data["time_series"], "pump")
            ts = data["time_series"]["pump"]
            pumps = values(filter(x -> x.first in keys(ts), data["pump"]))
            map(x -> x["flow_min"] = minimum(ts[string(x["index"])]["flow_min"]), pumps)
            map(x -> x["flow_min_forward"] = minimum(ts[string(x["index"])]["flow_min_forward"]), pumps)
            map(x -> x["flow_max"] = maximum(ts[string(x["index"])]["flow_max"]), pumps)
            map(x -> x["flow_max_reverse"] = maximum(ts[string(x["index"])]["flow_max_reverse"]), pumps)
            map(x -> x["z_min"] = minimum(ts[string(x["index"])]["z_min"]), pumps)
            map(x -> x["z_max"] = maximum(ts[string(x["index"])]["z_max"]), pumps)
        end
    end
end

function _relax_ne_pumps!(data::Dict{String,<:Any})
    if !_IM.ismultinetwork(data)
        if haskey(data, "time_series") && haskey(data["time_series"], "ne_pump")
            ts = data["time_series"]["ne_pump"]
            pumps = values(filter(x -> x.first in keys(ts), data["ne_pump"]))
            map(x -> x["flow_min"] = minimum(ts[string(x["index"])]["flow_min"]), pumps)
            map(x -> x["flow_min_forward"] = minimum(ts[string(x["index"])]["flow_min_forward"]), pumps)
            map(x -> x["flow_max"] = maximum(ts[string(x["index"])]["flow_max"]), pumps)
            map(x -> x["flow_max_reverse"] = maximum(ts[string(x["index"])]["flow_max_reverse"]), pumps)
            map(x -> x["z_min"] = minimum(ts[string(x["index"])]["z_min"]), pumps)
            map(x -> x["z_max"] = maximum(ts[string(x["index"])]["z_max"]), pumps)
        end
    end
end


function _set_pump_warm_start!(data::Dict{String, <:Any})
    for pump in values(data["pump"])
        flow_mid = 0.5 * (pump["flow_min"] + pump["flow_max"])

        pump["q_start"] = get(pump, "q", flow_mid)
        pump["qp_start"] = max(0.0, get(pump, "q", flow_mid))
        pump["qn_start"] = max(0.0, -get(pump, "q", flow_mid))

        pump["g_pump_start"] = get(pump, "g", 0.0)
        pump["P_pump_start"] = get(pump, "P", 0.0)

        pump["y_pump_start"] = get(pump, "q", 0.0) > 0.0 ? 1.0 : 0.0
        pump["z_pump_start"] = get(pump, "q", 0.0) > 0.0 ? 1.0 : 0.0
    end
end

function _set_ne_pump_warm_start!(data::Dict{String, <:Any})
    for pump in values(data["ne_pump"])
        flow_mid = 0.5 * (pump["flow_min"] + pump["flow_max"])

        pump["q_start"] = get(pump, "q", flow_mid)
        pump["qp_start"] = max(0.0, get(pump, "q", flow_mid))
        pump["qn_start"] = max(0.0, -get(pump, "q", flow_mid))

        pump["g_pump_start"] = get(pump, "g", 0.0)
        pump["P_pump_start"] = get(pump, "P", 0.0)

        pump["y_pump_start"] = get(pump, "q", 0.0) > 0.0 ? 1.0 : 0.0
        pump["z_pump_start"] = get(pump, "q", 0.0) > 0.0 ? 1.0 : 0.0
    end
end
