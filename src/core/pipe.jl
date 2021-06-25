function aggregate_pipes(subnetworks::Array{Dict{String, Any}, 1})
    return _aggregate_pipes(subnetworks, "pipe")
end


function aggregate_des_pipes(subnetworks::Array{Dict{String, Any}, 1})
    return _aggregate_pipes(subnetworks, "des_pipe")
end


function _aggregate_pipes(subnetworks::Array{Dict{String, Any}, 1}, pipe_type::String)
    pipes = deepcopy(subnetworks[1][pipe_type])

    for (i, x) in pipes
        x["flow_min"] = sum_subnetwork_values(subnetworks, pipe_type, i, "flow_min")
        x["flow_max"] = sum_subnetwork_values(subnetworks, pipe_type, i, "flow_max")
        x["flow_min_forward"] = sum_subnetwork_values(subnetworks, pipe_type, i, "flow_min_forward")
        x["flow_max_reverse"] = sum_subnetwork_values(subnetworks, pipe_type, i, "flow_max_reverse")
        x["status"] = any_subnetwork_values(subnetworks, pipe_type, i, "status")
    end

    return pipes
end


function correct_pipes!(data::Dict{String, <:Any})
    wm_data = get_wm_data(data)
    head_loss, viscosity = wm_data["head_loss"], wm_data["viscosity"]
    func! = x -> _correct_pipes!(x, head_loss, viscosity)
    apply_wm!(func!, data; apply_to_subnetworks = true)
end


function _correct_pipes!(data::Dict{String, <:Any}, head_loss::String, viscosity::Float64)
    capacity = _calc_capacity_max(data)
    base_length = get(data, "base_length", 1.0)
    base_time = get(data, "base_time", 1.0)

    for pipe in values(data["pipe"])
         # Get common connecting node data for later use.
        node_fr = data["node"][string(pipe["node_fr"])]
        node_to = data["node"][string(pipe["node_to"])]

        # Correct various pipe properties. The sequence is important, here.
        _correct_status!(pipe)
        _correct_flow_direction!(pipe)
        _correct_pipe_flow_bounds!(
            pipe, node_fr, node_to, head_loss,
            viscosity, capacity, base_length, base_time)
    end
end


function get_pipe_flow_lower_breakpoints(pipe::Dict{String, <:Any})
    if haskey(pipe, "flow_lower_breakpoints")
        return pipe["flow_lower_breakpoints"]
    else
        flow_min, flow_max = pipe["flow_min"], pipe["flow_max"]
        return [flow_min, flow_max]
    end
end


function get_pipe_flow_lower_breakpoints_positive(pipe::Dict{String, <:Any})
    lower_breakpoints = get_pipe_flow_lower_breakpoints(pipe)
    flows = filter(x -> x > 0.0, lower_breakpoints)
    lower_bound = max(0.0, get(pipe, "flow_min_forward", 0.0))
    flow_max = length(flows) > 0 ? maximum(flows) : lower_bound
    return lower_bound != flow_max ? vcat(lower_bound, flows) : [lower_bound]
end


function get_pipe_head_difference_lower_breakpoints_positive(
    pipe::Dict{String, <:Any}, head_loss::String, viscosity::Float64,
    base_length::Float64, base_time::Float64)
    flow_breakpoints = get_pipe_flow_lower_breakpoints_positive(pipe)
    exponent = _get_exponent_from_head_loss_form(head_loss)
    resistance = _calc_pipe_resistance(pipe, head_loss, viscosity, base_length, base_time)
    return (pipe["length"] * resistance) .* abs.(flow_breakpoints).^exponent
end


function get_pipe_flow_lower_breakpoints_negative(pipe::Dict{String, <:Any})
    lower_breakpoints = get_pipe_flow_lower_breakpoints(pipe)
    flows = filter(x -> x < 0.0, lower_breakpoints)
    upper_bound = min(0.0, get(pipe, "flow_max_reverse", 0.0))
    flow_min = length(flows) > 0 ? minimum(flows) : upper_bound
    return upper_bound != flow_min ? vcat(flows, upper_bound) : [upper_bound]
 end


function get_pipe_head_difference_lower_breakpoints_negative(
    pipe::Dict{String, <:Any}, head_loss::String, viscosity::Float64,
    base_length::Float64, base_time::Float64)
    flow_breakpoints = get_pipe_flow_lower_breakpoints_negative(pipe)
    exponent = _get_exponent_from_head_loss_form(head_loss)
    resistance = _calc_pipe_resistance(pipe, head_loss, viscosity, base_length, base_time)
    return (pipe["length"] * resistance) .* abs.(flow_breakpoints).^exponent
end


function get_pipe_flow_upper_breakpoints(pipe::Dict{String, <:Any})
    if haskey(pipe, "flow_upper_breakpoints")
        return pipe["flow_upper_breakpoints"]
    else
        flow_min, flow_max = pipe["flow_min"], pipe["flow_max"]
        return [flow_min, flow_max]
    end
end


function get_pipe_flow_upper_breakpoints_positive(pipe::Dict{String, <:Any})
    upper_breakpoints = get_pipe_flow_upper_breakpoints(pipe)
    flows = filter(x -> x > 0.0, upper_breakpoints)
    lower_bound = max(0.0, get(pipe, "flow_min_forward", 0.0))
    flow_max = length(flows) > 0 ? maximum(flows) : lower_bound
    return lower_bound != flow_max ? sort(vcat(lower_bound, flows)) : [lower_bound]
end


function get_pipe_head_difference_upper_breakpoints_positive(
    pipe::Dict{String, <:Any}, head_loss::String, viscosity::Float64,
    base_length::Float64, base_time::Float64)
    flow_breakpoints = get_pipe_flow_upper_breakpoints_positive(pipe)
    exponent = _get_exponent_from_head_loss_form(head_loss)
    resistance = _calc_pipe_resistance(pipe, head_loss, viscosity, base_length, base_time)
    return (pipe["length"] * resistance) .* abs.(flow_breakpoints).^exponent
end


function get_pipe_flow_upper_breakpoints_negative(pipe::Dict{String, <:Any})
    upper_breakpoints = get_pipe_flow_upper_breakpoints(pipe)
    flows = filter(x -> x < 0.0, upper_breakpoints)
    upper_bound = min(0.0, get(pipe, "flow_max_reverse", 0.0))
    flow_min = length(flows) > 0 ? minimum(flows) : upper_bound
    return upper_bound != flow_min ? sort(vcat(flows, upper_bound)) : [upper_bound]
end


function get_pipe_head_difference_upper_breakpoints_negative(
    pipe::Dict{String, <:Any}, head_loss::String, viscosity::Float64,
    base_length::Float64, base_time::Float64)
    flow_breakpoints = get_pipe_flow_upper_breakpoints_negative(pipe)
    exponent = _get_exponent_from_head_loss_form(head_loss)
    resistance = _calc_pipe_resistance(pipe, head_loss, viscosity, base_length, base_time)
    return (pipe["length"] * resistance) .* abs.(flow_breakpoints).^exponent
end


function correct_des_pipes!(data::Dict{String, <:Any})
    wm_data = get_wm_data(data)
    head_loss, viscosity = wm_data["head_loss"], wm_data["viscosity"]
    func! = x -> _correct_des_pipes!(x, head_loss, viscosity)
    apply_wm!(func!, data; apply_to_subnetworks = true)
end


function _correct_des_pipes!(data::Dict{String, <:Any}, head_loss::String, viscosity::Float64)
    capacity = _calc_capacity_max(data)
    base_length = get(data, "base_length", 1.0)
    base_time = get(data, "base_time", 1.0)

    for des_pipe in values(data["des_pipe"])
        # Get common connecting node data for later use.
        node_fr = data["node"][string(des_pipe["node_fr"])]
        node_to = data["node"][string(des_pipe["node_to"])]

        # Correct various pipe properties. The sequence is important, here.
        _correct_flow_direction!(des_pipe)
        _correct_pipe_flow_bounds!(des_pipe, node_fr, node_to,
            head_loss, viscosity, capacity, base_length, base_time)
    end
end


function set_pipe_breakpoints!(
    pipe::Dict{String, <:Any}, head_loss::String, viscosity::Float64, base_length::Float64,
    base_time::Float64, error_tolerance::Float64, length_tolerance::Float64)
    # Compute the product of pipe length and resistance.
    L_x_r = pipe["length"] * _calc_pipe_resistance(
        pipe, head_loss, viscosity, base_length, base_time)

    # Compute head loss function and derivative.
    exponent = uppercase(head_loss) == "H-W" ? 1.852 : 2.0
    f = x -> L_x_r * sign(x) * abs(x)^exponent
    f_dash = x -> exponent * L_x_r * (x * x)^(0.5 * exponent - 0.5)

    # Initialize the partitioning of flow breakpoints.
    partition = Vector{Float64}([pipe["flow_min"], pipe["flow_max"]])

    # Use PolyhedralRelaxations to determine partitions with desired accuracy.
    uvf_data = PolyhedralRelaxations.UnivariateFunctionData(
        f, f_dash, partition, error_tolerance, length_tolerance, 1.0e-6, 9e9, 0)
    PolyhedralRelaxations._refine_partition!(uvf_data)

    # Set the pipe lower and upper breakpoints using the above.
    pipe["flow_lower_breakpoints"] = uvf_data.partition
    pipe["flow_upper_breakpoints"] = uvf_data.partition
end


function _correct_pipe_flow_bounds!(
    pipe::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any},
    form::String, viscosity::Float64, capacity::Float64, base_length::Float64, base_time::Float64)
    # Calculate flow bounds, which depend on the head loss formulation.
    flow_min = _calc_pipe_flow_min(pipe, node_fr, node_to, form, viscosity, capacity, base_length, base_time)
    flow_max = _calc_pipe_flow_max(pipe, node_fr, node_to, form, viscosity, capacity, base_length, base_time)
    pipe["flow_min"], pipe["flow_max"] = flow_min, flow_max
    pipe["flow_min_forward"] = max(flow_min, get(pipe, "flow_min_forward", 0.0))
    pipe["flow_max_reverse"] = min(flow_max, get(pipe, "flow_max_reverse", 0.0))
end


function _calc_pipe_flow_min(
    pipe::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any},
    form::String, viscosity::Float64, capacity::Float64, base_length::Float64, base_time::Float64)
    # Calculate minimum flow bound depending on the head loss formulation.
    if uppercase(form) == "H-W"
        return _calc_pipe_flow_min_hw(pipe, node_fr, node_to, capacity, base_length, base_time)
    elseif uppercase(form) == "D-W"
        return _calc_pipe_flow_min_dw(pipe, node_fr, node_to, viscosity, capacity)
    end
end


function _calc_pipe_flow_max(
    pipe::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any},
    form::String, viscosity::Float64, capacity::Float64, base_length::Float64, base_time::Float64)
    # Calculate maximum flow bound depending on the head loss formulation.
    if uppercase(form) == "H-W"
        return _calc_pipe_flow_max_hw(pipe, node_fr, node_to, capacity, base_length, base_time)
    elseif uppercase(form) == "D-W"
        return _calc_pipe_flow_max_dw(pipe, node_fr, node_to, viscosity, capacity)
    end
end


function _calc_pipe_resistance(pipe::Dict{String, <:Any}, form::String, viscosity::Float64, base_length::Float64, base_time::Float64)
    if uppercase(form) == "H-W"
        return _calc_pipe_resistance_hw(pipe["diameter"], pipe["roughness"], base_length, base_time)
    elseif uppercase(form) == "D-W"
        diameter, roughness, speed = pipe["diameter"], pipe["roughness"], 10.0
        return _calc_pipe_resistance_dw(diameter, roughness, viscosity, speed, _DENSITY)
    end
end


function _calc_pipe_flow_min_hw(
    pipe::Dict{String, <:Any}, node_fr::Dict{String, <:Any},
    node_to::Dict{String, <:Any}, capacity::Float64, base_length::Float64, base_time::Float64)
    # Calculate minimum flow bound per the Hazen-Williams head loss formulation.
    resistance = _calc_pipe_resistance_hw(pipe["diameter"], pipe["roughness"], base_length, base_time)
    loss = get(node_fr, "head_min", -Inf) - get(node_to, "head_max", Inf)
    flow_min_loss = sign(loss) * (abs(loss) * inv(pipe["length"] * resistance))^inv(1.852)
    flow_min_dir = pipe["flow_direction"] == FLOW_DIRECTION_POSITIVE ? 0.0 : -Inf
    return max(-capacity, flow_min_loss, flow_min_dir, get(pipe, "flow_min", -Inf))
end


function _calc_pipe_flow_min_dw(
    pipe::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any},
    viscosity::Float64, capacity::Float64)
    # Calculate minimum flow bound per the Hazen-Williams head loss formulation.
    diameter, roughness, speed = pipe["diameter"], pipe["roughness"], 10.0
    resistance = _calc_pipe_resistance_dw(diameter, roughness, viscosity, speed, _DENSITY)
    loss = get(node_fr, "head_min", -Inf) - get(node_to, "head_max", Inf)
    flow_min_loss = sign(loss) * (abs(loss) * inv(pipe["length"] * resistance))^inv(2.0)
    flow_min_dir = pipe["flow_direction"] == FLOW_DIRECTION_POSITIVE ? 0.0 : -Inf
    return max(-capacity, flow_min_loss, flow_min_dir, get(pipe, "flow_min", -Inf))
end


function _calc_pipe_flow_max_hw(
    pipe::Dict{String, <:Any}, node_fr::Dict{String, <:Any},
    node_to::Dict{String, <:Any}, capacity::Float64, base_length::Float64, base_time::Float64)
    # Calculate maximum flow bound per the Hazen-Williams head loss formulation.
    resistance = _calc_pipe_resistance_hw(pipe["diameter"], pipe["roughness"], base_length, base_time)
    loss = get(node_fr, "head_max", Inf) - get(node_to, "head_min", -Inf)
    flow_max_loss = sign(loss) * (abs(loss) * inv(pipe["length"] * resistance))^inv(1.852)
    flow_max_dir = pipe["flow_direction"] == FLOW_DIRECTION_NEGATIVE ? 0.0 : Inf
    return min(capacity, flow_max_loss, flow_max_dir, get(pipe, "flow_max", Inf))
end


function _calc_pipe_flow_max_dw(
    pipe::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any},
    viscosity::Float64, capacity::Float64)
    # Calculate maximum flow bound per the Hazen-Williams head loss formulation.
    diameter, roughness, speed = pipe["diameter"], pipe["roughness"], 10.0
    resistance = _calc_pipe_resistance_dw(diameter, roughness, viscosity, speed, _DENSITY)
    loss = get(node_fr, "head_max", Inf) - get(node_to, "head_min", -Inf)
    flow_max_loss = sign(loss) * (abs(loss) * inv(pipe["length"] * resistance))^inv(2.0)
    flow_max_dir = pipe["flow_direction"] == FLOW_DIRECTION_NEGATIVE ? 0.0 : Inf
    return min(capacity, flow_max_loss, flow_max_dir, get(pipe, "flow_max", Inf))
end


function _get_exponent_from_head_loss_form(head_loss_form::String)
    return uppercase(head_loss_form) == "H-W" ? 1.852 : 2.0
end


"""
Computes the resistance for a pipe governed by the Hazen-Williams relationship.
"""
function _calc_pipe_resistance_hw(diameter::Float64, roughness::Float64, base_length::Float64, base_time::Float64)
    # Conversion factor for SI units. This value is in units of:
    # ((meters^3 / s)^1.852 / (meters^4.8704))^(1 / 1.852).
    k_si = 0.849 # This is the Hazen-Williams conversion factor given on Wikipedia.

    # Convert the conversion factor above (in SI units) to per-unit units.
    k = k_si * (((inv(base_length))^3)^1.852 / (inv(base_length))^4.8704)^(1 / 1.852) / inv(base_time)

    # Return the pipe resistance in the per-unit unit system.
    return 7.8828 * inv(k^1.852 * roughness^1.852 * diameter^4.8704)
end


"""
Computes the resistance for a pipe governed by the Darcy-Weisbach relationship.
"""
function _calc_pipe_resistance_dw(diameter::Float64, roughness::Float64, viscosity::Float64, speed::Float64, density::Float64)
    # Compute the Reynold's number of the fluid.
    reynolds_number = density * speed * diameter * inv(viscosity)

    # Use the same Colebrook approximation as in EPANET.
    w = 0.25 * pi * reynolds_number
    y_1 = 4.61841319859 * inv(w^0.9)
    y_2 = (roughness * inv(diameter)) * inv(3.7 * diameter) + y_1
    y_3 = -8.685889638e-01 * log(y_2)
    return 0.0826 * inv(diameter^5) * inv(y_3*y_3)
end


function _calc_head_loss_values(points::Array{Float64}, alpha::Float64)
    return [sign(x) * abs(x)^alpha for x in points]
end


function _calc_head_loss_oa(q::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, q_hat::Float64, exponent::Float64)
    return q_hat^exponent * z + exponent * q_hat^(exponent - 1.0) * (q - q_hat * z)
end


function set_pipe_warm_start!(data::Dict{String, <:Any})
    InfrastructureModels.apply!(_set_pipe_warm_start!, data, wm_it_name)
end


function _set_pipe_warm_start!(data::Dict{String, <:Any})
    for pipe in values(data["pipe"])
        flow_mid = 0.5 * (pipe["flow_min"] + pipe["flow_max"])
        
        pipe["q_start"] = get(pipe, "q", flow_mid)
        pipe["qp_start"] = max(0.0, get(pipe, "q", flow_mid))
        pipe["qn_start"] = max(0.0, -get(pipe, "q", flow_mid))
        pipe["y_pipe_start"] = get(pipe, "q", flow_mid) >= 0.0 ? 1.0 : 0.0

        node_fr = data["node"][string(pipe["node_fr"])]
        h_fr = get(node_fr, "h", 0.5 * (node_fr["head_min"] + node_fr["head_max"]))

        node_to = data["node"][string(pipe["node_to"])]
        h_to = get(node_to, "h", 0.5 * (node_to["head_min"] + node_to["head_max"]))

        pipe["dhp_start"] = max(0.0, h_fr - h_to)
        pipe["dhn_start"] = max(0.0, h_to - h_fr)
    end
end