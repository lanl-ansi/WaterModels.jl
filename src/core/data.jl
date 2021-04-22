# Functions for working with WaterModels data elements.


"Transform length values in SI units to per-unit units."
function _calc_length_per_unit_transform(data::Dict{String,<:Any})
    median_midpoint = _calc_node_head_median_midpoint(data)

    # Translation: convert number of meters to WaterModels-length.
    return x -> x / median_midpoint
end


"Correct flow direction attribute of edge-type components."
function _correct_flow_direction!(comp::Dict{String, <:Any})
    if !isa(comp["flow_direction"], FLOW_DIRECTION)
        comp["flow_direction"] = FLOW_DIRECTION(comp["flow_direction"])
    end
end


"Transform length values SI units to per-unit units."
function _calc_head_per_unit_transform(data::Dict{String,<:Any})
    median_midpoint = _calc_node_head_median_midpoint(data)

    # Translation: convert number of meters to WaterModels-length.
    return x -> x / median_midpoint
end


"Transform head values in per-unit units to SI units."
function _calc_head_per_unit_untransform(data::Dict{String,<:Any})
    median_midpoint = _calc_node_head_median_midpoint(data)
    return x -> x * median_midpoint + median_midpoint
end


"Transform time values in SI units to per-unit units."
function _calc_time_per_unit_transform(data::Dict{String,<:Any})
    # Translation: number of meters per WaterModels-length
    head_midpoint = _calc_node_head_median_midpoint(data)
    flow_midpoint = _calc_median_abs_flow_midpoint(data)

    # Translation: convert number of seconds to WaterModels-time.
    return x -> x / (head_midpoint^3 / flow_midpoint)
end

"Transform flow values in SI units to per-unit units."
function _calc_flow_per_unit_transform(data::Dict{String,<:Any})
    length_transform = _calc_length_per_unit_transform(data)
    time_transform = _calc_time_per_unit_transform(data)
    wm_vol_per_cubic_meter = length_transform(1.0)^3
    wm_time_per_second = time_transform(1.0)

    # Translation: convert cubic meters per second to WaterModels-flow.
    return x -> x * (wm_vol_per_cubic_meter / wm_time_per_second)
end


"Transform mass values in SI units to per-unit units."
function _calc_mass_per_unit_transform(data::Dict{String,<:Any})
    # Translation: convert kilograms to WaterModels-mass.
    return x -> x
end


function _calc_median_abs_flow_midpoint(data::Dict{String,<:Any})
    # Compute flow midpoints for all possible node-connecting components.
    q_des_pipe = [_calc_abs_flow_midpoint(x) for (i, x) in data["des_pipe"]]
    q_pipe = [_calc_abs_flow_midpoint(x) for (i, x) in data["pipe"]]
    q_pump = [_calc_abs_flow_midpoint(x) for (i, x) in data["pump"]]
    q_regulator = [_calc_abs_flow_midpoint(x) for (i, x) in data["regulator"]]
    q_short_pipe = [_calc_abs_flow_midpoint(x) for (i, x) in data["short_pipe"]]
    q_valve = [_calc_abs_flow_midpoint(x) for (i, x) in data["valve"]]
    q = vcat(q_des_pipe, q_pipe, q_pump, q_regulator, q_short_pipe, q_valve)

    # Return the median of all values computed above.
    return Statistics.median(q)
end


function _calc_abs_flow_midpoint(comp::Dict{String,Any})
    if comp["flow_direction"] in [0, UNKNOWN]
        qp_ub_mid = 0.5 * max(0.0, comp["flow_max"])
        qn_ub_mid = 0.5 * max(0.0, -comp["flow_min"])
        return max(qp_ub_mid, qn_ub_mid)
    elseif comp["flow_direction"] in [1, POSITIVE]
        return 0.5 * max(0.0, comp["flow_max"])
    elseif comp["flow_direction"] in [-1, NEGATIVE]
        return 0.5 * max(0.0, -comp["flow_min"])
    end
end


"""
Turns given single network data into multinetwork data with `count` replicates of the given
network. Note that this function performs a deepcopy of the network data. Significant
multinetwork space savings can often be achieved by building application specific methods of
building a multinetwork with minimal data replication (e.g., through storing references).
"""
function replicate(data::Dict{String,<:Any}, count::Int; global_keys::Set{String}=Set{String}())
    return _IM.replicate(data, count, union(global_keys, _wm_global_keys))
end


function _remove_last_networks!(data::Dict{String, <:Any}; last_num_steps::Int = length(nw_ids(wm)))
    @assert _IM.ismultinetwork(data) # Ensure the data being operated on is multinetwork.
    network_ids = reverse(sort([parse(Int, x) for x in keys(data["nw"])]))
    network_ids_exclude = network_ids[1:min(length(network_ids), last_num_steps)]
    network_ids_exclude_str = [string(x) for x in network_ids_exclude]
    data["nw"] = filter(x -> !(x.first in network_ids_exclude_str), data["nw"])
end


function _calc_head_max(data::Dict{String, <:Any})
    # Compute the maximum elevation of all nodes in the network.
    head_max = maximum(get(node, "head_min", -Inf) for (i, node) in data["node"])
    head_max = maximum(get(node, "head_nominal", head_max) for (i, node) in data["node"])
    head_max = maximum(get(node, "head_max", head_max) for (i, node) in data["node"])

    for (i, tank) in data["tank"]
        # Consider maximum tank head in computation of head_max.
        elevation = data["node"][string(tank["node"])]["elevation"]
        head_max = max(head_max, elevation + tank["max_level"])
    end

    for (i, pump) in data["pump"]
        # Consider possible pump head gains in computation of head_max.
        node_fr = data["node"][string(pump["node_fr"])]
        node_to = data["node"][string(pump["node_to"])]
        head_gain = _calc_pump_head_gain_max(pump, node_fr, node_to)
        head_max = max(head_max, node_to["elevation"] + head_gain)
    end

    for (i, regulator) in data["regulator"]
        # Consider possible downstream regulator heads in computation of head_max.
        p_setting, node_to_index = regulator["setting"], string(regulator["node_to"])
        elevation = data["node"][node_to_index]["elevation"]
        head_max = max(head_max, elevation + p_setting)
    end

    return head_max
end


function _calc_capacity_max(data::Dict{String, <:Any})
    # Include the sum of all maximal flows from demands.
    capacity = sum(x["flow_max"] for (i, x) in data["demand"])

    for (i, tank) in data["tank"]
        # Add the sum of maximum possible demanded flow from tanks.
        surface_area = 0.25 * pi * tank["diameter"]^2
        volume_min = max(tank["min_vol"], surface_area * tank["min_level"])
        volume_max = surface_area * tank["max_level"]
        capacity += (volume_max - volume_min) * inv(data["time_step"])
    end

    # Return the maximum capacity of the network.
    return capacity
end


"Turns a single network with a `time_series` data block into a multinetwork."
function make_multinetwork(data::Dict{String, <:Any}; global_keys::Set{String} = Set{String}())
    return _IM.make_multinetwork(data, wm_it_name, union(global_keys, _wm_global_keys))
end


function split_multinetwork(data::Dict{String, <:Any}, nw_ids::Array{Array{String, 1}, 1})
    # Ensure the data has the multinetwork attribute.
    @assert _IM.ismultinetwork(data) == true

    # Get sub-multinetwork datasets indexed by the "nw" key.
    sub_mn = [Dict{String, Any}(i => deepcopy(data["nw"][i]) for i in ids) for ids in nw_ids]

    # Get all data not associated with multinetwork components
    g_data = filter(x -> x.first != "nw", data)

    # Return the new, split multinetwork data dictionaries.
    return [merge(deepcopy(g_data), Dict{String, Any}("nw" => sub_mn[i])) for i in 1:length(nw_ids)]
end


function set_start!(data::Dict{String,<:Any}, component_type::String, var_name::String, start_name::String)
    if _IM.ismultinetwork(data)
        for (n, nw) in data["nw"]
            comps = values(nw[component_type])
            map(x -> x[start_name] = x[var_name], comps)
        end
    else
        comps = values(data[component_type])
        map(x -> x[start_name] = x[var_name], comps)
    end
end


function set_flow_start!(data::Dict{String,<:Any})
    set_start!(data, "pipe", "q", "q_pipe_start")
    set_start!(data, "pump", "q", "q_pump_start")
    set_start!(data, "regulator", "q", "q_regulator_start")
    set_start!(data, "short_pipe", "q", "q_short_pipe_start")
    set_start!(data, "valve", "q", "q_valve_start")
    set_start!(data, "reservoir", "q", "q_reservoir_start")
    set_start!(data, "tank", "q", "q_tank_start")
end


function set_head_start!(data::Dict{String,<:Any})
    set_start!(data, "node", "h", "h_start")
    set_start!(data, "pump", "g", "g_pump_start")
end


function set_indicator_start!(data::Dict{String,<:Any})
    set_start!(data, "pump", "status", "z_pump_start")
    set_start!(data, "regulator", "status", "z_regulator_start")
    set_start!(data, "valve", "status", "z_valve_start")
end


function set_start_all!(data::Dict{String,<:Any})
    set_flow_start!(data)
    set_head_start!(data)
    set_indicator_start!(data)
end


function fix_all_flow_directions!(data::Dict{String,<:Any})
    fix_flow_directions!(data, "pipe")
    fix_flow_directions!(data, "short_pipe")
    fix_flow_directions!(data, "pump")
    fix_flow_directions!(data, "regulator")
    fix_flow_directions!(data, "valve")
end


function fix_flow_directions!(data::Dict{String,<:Any}, component_type::String)
    if _IM.ismultinetwork(data)
        for (n, nw) in data["nw"]
            comps = values(nw[component_type])
            _fix_flow_direction!.(comps)
        end
    else
        comps = values(data[component_type])
        _fix_flow_direction!.(comps)
    end
end


function _fix_flow_direction!(component::Dict{String,<:Any})
    if haskey(component, "q") && !isapprox(component["q"], 0.0; atol=1.0e-6)
        component["y_min"] = component["q"] > 0.0 ? 1.0 : 0.0
        component["y_max"] = component["q"] > 0.0 ? 1.0 : 0.0
    end
end


function fix_all_indicators!(data::Dict{String,<:Any})
    fix_indicators!(data, "pump")
    fix_indicators!(data, "regulator")
    fix_indicators!(data, "valve")
end

function fix_indicators!(data::Dict{String,<:Any}, component_type::String)
    if _IM.ismultinetwork(data)
        for (n, nw) in data["nw"]
            comps = values(nw[component_type])
            _fix_indicator!.(comps)
        end
    else
        comps = values(data[component_type])
        _fix_indicator!.(comps)
    end
end


function _fix_indicator!(component::Dict{String,<:Any})
    if haskey(component, "status") # Assumes this is binary.
        component["z_min"] = component["status"] > 0.5 ? 1.0 : 0.0
        component["z_max"] = component["status"] > 0.5 ? 1.0 : 0.0
    end
end


"Translate a multinetwork dataset to a snapshot dataset with dispatchable components."
function _relax_network!(data::Dict{String,<:Any})
    _relax_nodes!(data)
    _relax_tanks!(data)
    _relax_reservoirs!(data)
    _relax_demands!(data)
end


"Convenience function for recomputing component bounds, e.g., after modifying data."
function recompute_bounds!(data::Dict{String, <:Any})
    # Clear the existing flow bounds for node-connecting components.
    for comp_type in ["pipe", "des_pipe", "pump", "regulator", "short_pipe", "valve"]
        map(x -> x["flow_min"] = -Inf, values(data[comp_type]))
        map(x -> x["flow_max"] = Inf, values(data[comp_type]))
    end

    # Recompute bounds and correct data.
    correct_network_data!(data)
end


function sum_subnetwork_values(subnetworks::Array{Dict{String, Any},1}, comp_name::String, index::String, key::String)
    return sum(nw[comp_name][index][key] for nw in subnetworks)
end


function max_subnetwork_values(subnetworks::Array{Dict{String, Any},1}, comp_name::String, index::String, key::String)
    return maximum(nw[comp_name][index][key] for nw in subnetworks)
end


function min_subnetwork_values(subnetworks::Array{Dict{String, Any},1}, comp_name::String, index::String, key::String)
    return minimum(nw[comp_name][index][key] for nw in subnetworks)
end


function all_subnetwork_values(subnetworks::Array{Dict{String, Any},1}, comp_name::String, index::String, key::String)
    return all(nw[comp_name][index][key] for nw in subnetworks)
end


function any_subnetwork_values(subnetworks::Array{Dict{String, Any},1}, comp_name::String, index::String, key::String)
    return any(nw[comp_name][index][key] == 1 for nw in subnetworks) ? 1 : 0
end


function aggregate_time_steps(subnetworks::Array{Dict{String, Any}, 1})
    return sum(x["time_step"] for x in subnetworks)
end


function aggregate_subnetworks(data::Dict{String, Any}, nw_ids::Array{String, 1})
    # Initialize important metadata and dictionary.
    subnetworks = [data["nw"][x] for x in nw_ids]
    time_step = aggregate_time_steps(subnetworks)
    data_agg = Dict{String,Any}("time_step" => time_step)

    # Aggregate nodal components.
    data_agg["node"] = aggregate_nodes(subnetworks)
    data_agg["demand"] = aggregate_demands(subnetworks)
    data_agg["reservoir"] = aggregate_reservoirs(subnetworks)
    data_agg["tank"] = aggregate_tanks(subnetworks)

    # Aggregate node-connecting componets.
    data_agg["pipe"] = aggregate_pipes(subnetworks)
    data_agg["des_pipe"] = aggregate_des_pipes(subnetworks)
    data_agg["pump"] = aggregate_pumps(subnetworks)
    data_agg["regulator"] = aggregate_regulators(subnetworks)
    data_agg["short_pipe"] = aggregate_short_pipes(subnetworks)
    data_agg["valve"] = aggregate_valves(subnetworks)

    # Return the time-aggregated network.
    return data_agg
end

function make_temporally_aggregated_multinetwork(data::Dict{String, <:Any}, nw_ids::Array{Array{String, 1}, 1})
    # Initialize the temporally aggregated multinetwork.
    new_nw_ids = [string(i) for i in 1:length(nw_ids)]
    data_agg = Dict{String, Any}("nw" => Dict{String, Any}())

    for n in 1:length(new_nw_ids)
        # Construct and add the aggregation of subnetworks.
        new_nw_id, old_nw_ids = new_nw_ids[n], nw_ids[n]
        data_agg["nw"][new_nw_id] = aggregate_subnetworks(data, old_nw_ids)
    end

    data_agg["name"], data_agg["per_unit"] = data["name"], data["per_unit"]
    data_agg["viscosity"], data_agg["multinetwork"] = data["viscosity"], true
    duration = sum(x["time_step"] for (i, x) in data_agg["nw"])
    data_agg["duration"], data_agg["head_loss"] = duration, data["head_loss"]

    # Return the temporally aggregated multinetwork.
    return data_agg
end


function _transform_flow_key!(comp::Dict{String,<:Any}, key::String, transform_flow::Function)
    haskey(comp, key) && (comp[key] = transform_flow(comp[key]))
end


function _transform_flows!(comp::Dict{String,<:Any}, transform_flow::Function)
    _transform_flow_key!(comp, "flow_min", transform_flow)
    _transform_flow_key!(comp, "flow_max", transform_flow)
    _transform_flow_key!(comp, "flow_min_forward", transform_flow)
    _transform_flow_key!(comp, "flow_max_reverse", transform_flow)
    _transform_flow_key!(comp, "minor_loss", transform_flow)
end


function _make_per_unit_nodes!(data::Dict{String,<:Any}, transform_head::Function)
    for (i, node) in data["node"]
        node["elevation"] = transform_head(node["elevation"])
        node["head_min"] = transform_head(node["head_min"])
        node["head_max"] = transform_head(node["head_max"])
        node["head_nominal"] = transform_head(node["head_nominal"])
    end
end


function _make_per_unit_reservoir!(data::Dict{String,<:Any}, transform_head::Function)
    for (i, reservoir) in data["reservoir"]
        reservoir["head_nominal"] = transform_head(reservoir["head_nominal"])
    end
end


function _make_per_unit_demands!(data::Dict{String,<:Any}, transform_flow::Function)
    for (i, demand) in data["demand"]
        demand["flow_min"] = transform_flow(demand["flow_min"])
        demand["flow_max"] = transform_flow(demand["flow_max"])
        demand["flow_nominal"] = transform_flow(demand["flow_nominal"])
    end
end


function _make_per_unit_tanks!(data::Dict{String,<:Any}, transform_length::Function)
    for (i, tank) in data["tank"]
        tank["min_level"] = transform_length(tank["min_level"])
        tank["max_level"] = transform_length(tank["max_level"])
        tank["init_level"] = transform_length(tank["init_level"])
        tank["diameter"] = transform_length(tank["diameter"])
        tank["min_vol"] = transform_length(tank["min_vol"])^3
    end
end


function _make_per_unit_heads!(data::Dict{String,<:Any}, transform_head::Function)
    _make_per_unit_nodes!(data, transform_head)
end


function _make_per_unit_flows!(data::Dict{String,<:Any}, transform_flow::Function)
    for type in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
        map(x -> _transform_flows!(x, transform_flow), values(data[type]))
    end
end


function _make_per_unit_pipes!(data::Dict{String,<:Any}, transform_length::Function)
    for (i, pipe) in data["pipe"]
        pipe["length"] = transform_length(pipe["length"])
        pipe["diameter"] = transform_length(pipe["diameter"])
    end
end


function _make_per_unit_des_pipes!(data::Dict{String,<:Any}, transform_length::Function)
    for (i, des_pipe) in data["des_pipe"]
        des_pipe["length"] = transform_length(des_pipe["length"])
        des_pipe["diameter"] = transform_length(des_pipe["diameter"])
    end
end


function _make_per_unit_pumps!(data::Dict{String,<:Any}, transform_flow::Function, transform_length::Function)
    for (i, pump) in data["pump"]
        pump["head_curve"] = [(transform_flow(x[1]), x[2]) for x in pump["head_curve"]]
        pump["head_curve"] = [(x[1], transform_length(x[2])) for x in pump["head_curve"]]

        if haskey(pump, "efficiency_curve")
            pump["efficiency_curve"] = [(transform_flow(x[1]), x[2]) for x in pump["efficiency_curve"]]
        end
    end
end


function _calc_scaled_gravity(base_length::Float64, base_time::Float64)
    return _GRAVITY * base_length / (base_time)^2
end


function make_per_unit!(data::Dict{String,<:Any})
    if data["per_unit"] == false
        # Precompute per-unit transformation functions.
        mass_transform = _calc_mass_per_unit_transform(data)
        length_transform = _calc_length_per_unit_transform(data)
        head_transform = _calc_head_per_unit_transform(data)
        flow_transform = _calc_flow_per_unit_transform(data)
        time_transform = _calc_time_per_unit_transform(data)

        # Apply per-unit transformations to node-connecting components.
        _make_per_unit_flows!(data, flow_transform)
        _make_per_unit_pipes!(data, length_transform)
        _make_per_unit_des_pipes!(data, length_transform)
        _make_per_unit_pumps!(data, flow_transform, length_transform)

        # Apply per-unit transformations to nodal components.
        _make_per_unit_nodes!(data, head_transform)
        _make_per_unit_demands!(data, flow_transform)
        _make_per_unit_tanks!(data, length_transform)
        _make_per_unit_reservoir!(data, head_transform)

        # Convert viscosity to per-unit units.
        data["viscosity"] *= mass_transform(1.0) /
            (length_transform(1.0) * time_transform(1.0))

        # Store important base values for later conversions.
        data["base_length"] = length_transform(1.0)
        data["base_time"] = time_transform(1.0)

        # Set the per-unit flag.
        data["time_step"] = time_transform(data["time_step"])
        data["per_unit"] = true
    end
end