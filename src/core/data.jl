# Functions for working with WaterModels data elements.

"WaterModels wrapper for the InfrastructureModels `apply!` function."
function apply_wm!(func!::Function, data::Dict{String, <:Any}; apply_to_subnetworks::Bool = true)
    _IM.apply!(func!, data, wm_it_name; apply_to_subnetworks = apply_to_subnetworks)
end


function correct_enums!(data::Dict{String,<:Any})
    correct_statuses!(data)
    correct_flow_directions!(data)
    correct_pump_head_curve_forms!(data)
end


function correct_flow_directions!(data::Dict{String,<:Any})
    apply_wm!(_correct_flow_directions!, data; apply_to_subnetworks = true)
end


function _correct_flow_directions!(data::Dict{String,<:Any})
    for component_type in ["pipe", "des_pipe", "short_pipe", "pump", "regulator", "valve"]
        components = values(get(data, component_type, Dict{String,Any}()))
        _correct_flow_direction!.(components)
    end
end


function correct_statuses!(data::Dict{String,<:Any})
    apply_wm!(_correct_statuses!, data; apply_to_subnetworks = true)
end


function _correct_statuses!(data::Dict{String,<:Any})
    edge_types = ["pipe", "des_pipe", "short_pipe", "pump", "regulator", "valve"]
    node_types = ["node", "demand", "reservoir", "tank", "reservoir"]

    for component_type in vcat(edge_types, node_types)
        components = values(get(data, component_type, Dict{String,Any}()))
        _correct_status!.(components)
    end
end


"WaterModels wrapper for the InfrastructureModels `get_data` function."
function get_data_wm(func::Function, data::Dict{String, <:Any}; apply_to_subnetworks::Bool = true)
    return _IM.get_data(func, data, wm_it_name; apply_to_subnetworks = apply_to_subnetworks)
end


"Convenience function for retrieving the water-only portion of network data."
function get_wm_data(data::Dict{String, <:Any})
    return _IM.ismultiinfrastructure(data) ? data["it"][wm_it_name] : data
end


function set_flow_partitions_num!(data::Dict{String, <:Any}, num_points::Int)
    wm_data = get_wm_data(data)
    func! = x -> _set_flow_partitions_num!(x, num_points)
    apply_wm!(func!, wm_data; apply_to_subnetworks = true)
end


function _set_flow_partitions_num!(data::Dict{String, <:Any}, num_points::Int)
    for pipe in values(get(data, "pipe", Dict{String, Any}()))
        flow_min, flow_max = pipe["flow_min"], pipe["flow_max"]

        if flow_min < flow_max
            partition = range(flow_min, flow_max; length = num_points)
            pipe["flow_partition"] = collect(partition)
        else
            pipe["flow_partition"] = [flow_min]
        end
    end

    for pump in values(get(data, "pump", Dict{String, Any}()))
        flow_min, flow_max = pump["flow_min_forward"], pump["flow_max"]

        if flow_min < flow_max
            partition = range(flow_min, flow_max; length = num_points)
            pump["flow_partition"] = collect(partition)
        else
            pump["flow_partition"] = [flow_min]
        end
    end
end


function set_flow_partitions_si!(data::Dict{String, <:Any}, error_tolerance::Float64, length_tolerance::Float64)
    wm_data = get_wm_data(data)
    head_loss, viscosity = wm_data["head_loss"], wm_data["viscosity"]
    base_flow = wm_data["per_unit"] ? wm_data["base_flow"] : 1.0
    base_length = wm_data["per_unit"] ? wm_data["base_length"] : 1.0
    base_mass = wm_data["per_unit"] ? wm_data["base_mass"] : 1.0
    base_time = wm_data["per_unit"] ? wm_data["base_time"] : 1.0

    func! = x -> _set_flow_partitions_si!(
        x, error_tolerance / base_length, length_tolerance / base_flow,
        head_loss, viscosity, base_length, base_mass, base_time)

    apply_wm!(func!, wm_data; apply_to_subnetworks = true)
end


function _set_flow_partitions_si!(
    data::Dict{String, <:Any}, error_tolerance::Float64, length_tolerance::Float64,
    head_loss::String, viscosity::Float64, base_length::Float64, base_mass::Float64, base_time::Float64)
    # Set partitions for all pipes in the network.
    for pipe in values(get(data, "pipe", Dict{String, Any}()))
        set_pipe_flow_partition!(
            pipe, head_loss, viscosity, base_length, base_mass,
            base_time, error_tolerance, length_tolerance)
    end

    # Set partitions for all design pipes in the network.
    for des_pipe in values(get(data, "des_pipe", Dict{String, Any}()))
        set_pipe_flow_partition!(
            des_pipe, head_loss, viscosity, base_length, base_mass,
            base_time, error_tolerance, length_tolerance)
    end

    # Set partitions for all pumps in the network.
    for pump in values(get(data, "pump", Dict{String, Any}()))
        set_pump_flow_partition!(pump, error_tolerance, length_tolerance)
    end
end


"Get function to transform length values in SI units to per-unit units."
function _calc_length_per_unit_transform(data::Dict{String,<:Any})
    wm_data = get_wm_data(data)

    if haskey(wm_data, "base_length")
        @assert wm_data["base_length"] > 0.0
        return x -> x / wm_data["base_length"]
    else
        median_midpoint = _calc_node_head_median_midpoint(wm_data)
        @assert median_midpoint > 0.0
        return x -> x / median_midpoint
    end
end


"Get function to transform length values in per-unit units to SI units."
function _calc_length_si_transform(data::Dict{String,<:Any})
    wm_data = get_wm_data(data)
    @assert haskey(wm_data, "base_length")
    return x -> x * wm_data["base_length"]
end


"Correct flow direction attribute of edge-type components."
function _correct_flow_direction!(comp::Dict{String, <:Any})
    flow_direction = get(comp, "flow_direction", FLOW_DIRECTION_UNKNOWN)

    if isa(flow_direction, FLOW_DIRECTION)
        comp["flow_direction"] = flow_direction
    else
        comp["flow_direction"] = FLOW_DIRECTION(flow_direction)
    end
end


"Correct status attribute of a component."
function _correct_status!(comp::Dict{String, <:Any})
    status = get(comp, "status", STATUS_UNKNOWN)

    if isa(status, STATUS)
        comp["status"] = status
    else
        comp["status"] = STATUS(Int(status))
    end
end


"Get function to transform length values SI units to per-unit units."
function _calc_head_per_unit_transform(data::Dict{String,<:Any})
    wm_data = get_wm_data(data)

    if haskey(wm_data, "base_head")
        @assert wm_data["base_head"] > 0.0
        return x -> x / wm_data["base_head"]
    else
        median_midpoint = _calc_node_head_median_midpoint(wm_data)
        @assert median_midpoint > 0.0
        return x -> x / median_midpoint
    end
end


"Get function to transform head values in per-unit units to SI units."
function _calc_head_si_transform(data::Dict{String,<:Any})
    wm_data = get_wm_data(data)
    @assert haskey(wm_data, "base_head")
    return x -> x * wm_data["base_head"]
end


"Get function to transform time values in SI units to per-unit units."
function _calc_time_per_unit_transform(data::Dict{String,<:Any})
    wm_data = get_wm_data(data)

    if haskey(wm_data, "base_time")
        @assert wm_data["base_time"] > 0.0
        return x -> x / wm_data["base_time"]
    else
        # Translation: number of meters per WaterModels length.
        head_midpoint = _calc_node_head_median_midpoint(wm_data)

        # Translation: number of cubic meters per second per WaterModels flow.
        flow_midpoint = _calc_median_abs_flow_midpoint(wm_data)

        # Translation: convert number of seconds to WaterModels time.
        return x -> x / (head_midpoint^3 / flow_midpoint)
    end
end


"Get function to transform time values in per-unit units to SI units."
function _calc_time_si_transform(data::Dict{String,<:Any})
    wm_data = get_wm_data(data)
    @assert haskey(wm_data, "base_time")
    return x -> x * wm_data["base_time"]
end

"Get function to transform flow values in SI units to per-unit units."
function _calc_flow_per_unit_transform(data::Dict{String,<:Any})
    wm_data = get_wm_data(data)

    if haskey(wm_data, "base_flow")
        return x -> x / wm_data["base_flow"]
    else
        length_transform = _calc_length_per_unit_transform(wm_data)
        time_transform = _calc_time_per_unit_transform(wm_data)
        wm_vol_per_cubic_meter = length_transform(1.0)^3
        wm_time_per_second = time_transform(1.0)

        # Translation: convert cubic meters per second to WaterModels-flow.
        return x -> x * (wm_vol_per_cubic_meter / wm_time_per_second)
    end
end


"Get function to transform flow values in per-unit units to SI units."
function _calc_flow_si_transform(data::Dict{String,<:Any})
    wm_data = get_wm_data(data)
    @assert haskey(wm_data, "base_flow")
    return x -> x * wm_data["base_flow"]
end


"Get function to transform mass values in SI units to per-unit units."
function _calc_mass_per_unit_transform(data::Dict{String,<:Any})
    # Translation: convert kilograms to WaterModels-mass.
    return x -> x / get(get_wm_data(data), "base_mass", 1.0e19)
end


"Get function to transform mass values in per-unit units to SI units."
function _calc_mass_si_transform(data::Dict{String,<:Any})
    wm_data = get_wm_data(data)
    @assert haskey(wm_data, "base_mass")
    return x -> x * wm_data["base_mass"]
end


function _calc_median_abs_flow_midpoint(data::Dict{String,<:Any})
    wm_data = get_wm_data(data)

    if ismultinetwork(data)
        flow_midpoints = Vector{Float64}([])

        for (n, nw) in data["nw"]
            q_des_pipe = [_calc_abs_flow_midpoint(x) for (i, x) in nw["des_pipe"]]
            q_pipe = [_calc_abs_flow_midpoint(x) for (i, x) in nw["pipe"]]
            q_pump = [_calc_abs_flow_midpoint(x) for (i, x) in nw["pump"]]
            q_regulator = [_calc_abs_flow_midpoint(x) for (i, x) in nw["regulator"]]
            q_short_pipe = [_calc_abs_flow_midpoint(x) for (i, x) in nw["short_pipe"]]
            q_valve = [_calc_abs_flow_midpoint(x) for (i, x) in nw["valve"]]
            q = vcat(q_des_pipe, q_pipe, q_pump, q_regulator, q_short_pipe, q_valve)
            append!(flow_midpoints, q)
        end
    else
        q_des_pipe = [_calc_abs_flow_midpoint(x) for (i, x) in wm_data["des_pipe"]]
        q_pipe = [_calc_abs_flow_midpoint(x) for (i, x) in wm_data["pipe"]]
        q_pump = [_calc_abs_flow_midpoint(x) for (i, x) in wm_data["pump"]]
        q_regulator = [_calc_abs_flow_midpoint(x) for (i, x) in wm_data["regulator"]]
        q_short_pipe = [_calc_abs_flow_midpoint(x) for (i, x) in wm_data["short_pipe"]]
        q_valve = [_calc_abs_flow_midpoint(x) for (i, x) in wm_data["valve"]]
        flow_midpoints = vcat(q_des_pipe, q_pipe, q_pump, q_regulator, q_short_pipe, q_valve)
    end

    # Return the median of all values computed above.
    return Statistics.median(flow_midpoints)
end


function _calc_abs_flow_midpoint(comp::Dict{String,Any})
    if comp["flow_direction"] in [0, FLOW_DIRECTION_UNKNOWN]
        qp_ub_mid = 0.5 * max(0.0, comp["flow_max"])
        qn_ub_mid = 0.5 * max(0.0, -comp["flow_min"])
        return max(qp_ub_mid, qn_ub_mid)
    elseif comp["flow_direction"] in [1, FLOW_DIRECTION_POSITIVE]
        return 0.5 * max(0.0, comp["flow_max"])
    elseif comp["flow_direction"] in [-1, FLOW_DIRECTION_NEGATIVE]
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


function _calc_head_offset(data::Dict{String, <:Any})
    wm_data = get_wm_data(data)

    if wm_data["per_unit"]
        head_transform = _calc_head_per_unit_transform(wm_data)
        return head_transform(-100.0) # -100 meters of head, scaled.
    else
        return -100.0 # -100 meters of head, unscaled.
    end
end


function _calc_head_max(data::Dict{String, <:Any})
    wm_data = get_wm_data(data)

    # Compute the maximum elevation of all nodes in the network.
    head_max = maximum(get(node, "head_min", -Inf) for (i, node) in wm_data["node"])
    head_max = maximum(get(node, "head_nominal", head_max) for (i, node) in wm_data["node"])
    head_max = maximum(get(node, "head_max", head_max) for (i, node) in wm_data["node"])

    for (i, tank) in wm_data["tank"]
        # Consider maximum tank head in computation of head_max.
        elevation = wm_data["node"][string(tank["node"])]["elevation"]
        head_max = max(head_max, elevation + tank["max_level"])
    end

    for (i, pump) in wm_data["pump"]
        # Consider possible pump head gains in computation of head_max.
        node_fr = wm_data["node"][string(pump["node_fr"])]
        node_to = wm_data["node"][string(pump["node_to"])]
        head_gain = _calc_pump_head_gain_max(pump, node_fr, node_to)
        head_max = max(head_max, node_to["elevation"] + head_gain)
    end

    for (i, regulator) in wm_data["regulator"]
        # Consider possible downstream regulator heads in computation of head_max.
        p_setting, node_to_index = regulator["setting"], string(regulator["node_to"])
        elevation = wm_data["node"][node_to_index]["elevation"]
        head_max = max(head_max, elevation + p_setting)
    end

    return head_max
end


function _calc_capacity_max(data::Dict{String, <:Any})
    wm_data = get_wm_data(data)

    # Include the sum of all maximal flows from demands.
    capacity = sum(x["flow_max"] for (i, x) in wm_data["demand"])

    for (i, tank) in wm_data["tank"]
        # Add the sum of maximum possible demanded flow from tanks.
        surface_area = 0.25 * pi * tank["diameter"]^2
        volume_min = max(tank["min_vol"], surface_area * tank["min_level"])
        volume_max = surface_area * tank["max_level"]
        capacity += (volume_max - volume_min) / wm_data["time_step"]
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
    sub_mn = [Dict{String, Any}(i => deepcopy(data["nw"][i])
        for i in ids) for ids in nw_ids]

    # Get all data not associated with multinetwork components
    g_data = filter(x -> x.first != "nw", data)

    # Return the new, split multinetwork data dictionaries.
    return [merge(deepcopy(g_data), Dict{String, Any}("nw" =>
        sub_mn[i])) for i in 1:length(nw_ids)]
end


function set_start!(data::Dict{String,<:Any}, component_type::String, var_name::String, start_name::String)
    wm_data = get_wm_data(data)

    if _IM.ismultinetwork(wm_data)
        for (n, nw) in wm_data["nw"]
            comps = values(nw[component_type])
            map(x -> x[start_name] = x[var_name], comps)
        end
    else
        comps = values(wm_data[component_type])
        map(x -> x[start_name] = x[var_name], comps)
    end
end


function set_direction_start_from_flow!(data::Dict{String,<:Any}, component_type::String, var_name::String, start_name::String)
    wm_data = get_wm_data(data)

    if _IM.ismultinetwork(wm_data)
        for (n, nw) in wm_data["nw"]
            comps = values(nw[component_type])
            map(x -> x[start_name] = x[var_name] > 0.0 ? 1 : 0, comps)
        end
    else
        comps = values(wm_data[component_type])
        map(x -> x[start_name] = x[var_name] > 0.0 ? 1 : 0, comps)
    end
end


function set_flow_start!(data::Dict{String,<:Any})
    apply_wm!(_set_flow_start!, data)
end


function _set_flow_start!(data::Dict{String,<:Any})
    set_start!(data, "pipe", "q", "q_pipe_start")
    set_start!(data, "pump", "q", "q_pump_start")
    set_start!(data, "regulator", "q", "q_regulator_start")
    set_start!(data, "short_pipe", "q", "q_short_pipe_start")
    set_start!(data, "valve", "q", "q_valve_start")
    set_start!(data, "reservoir", "q", "q_reservoir_start")
    set_start!(data, "tank", "q", "q_tank_start")
end


function set_flow_direction_start!(data::Dict{String,<:Any})
    apply_wm!(_set_flow_direction_start!, data)
end


function _set_flow_direction_start!(data::Dict{String,<:Any})
    set_direction_start_from_flow!(data, "pipe", "q", "y_pipe_start")
    set_direction_start_from_flow!(data, "pump", "q", "y_pump_start")
    set_direction_start_from_flow!(data, "regulator", "q", "y_regulator_start")
    set_direction_start_from_flow!(data, "short_pipe", "q", "y_short_pipe_start")
    set_direction_start_from_flow!(data, "valve", "q", "y_valve_start")
end


function set_head_start!(data::Dict{String,<:Any})
    apply_wm!(_set_head_start!, data)
end


function _set_head_start!(data::Dict{String,<:Any})
    set_start!(data, "node", "h", "h_start")
    set_start!(data, "pump", "g", "g_pump_start")
end


function set_indicator_start!(data::Dict{String,<:Any})
    apply_wm!(_set_indicator_start!, data)
end


function _set_indicator_start!(data::Dict{String,<:Any})
    set_start!(data, "pump", "status", "z_pump_start")
    set_start!(data, "regulator", "status", "z_regulator_start")
    set_start!(data, "valve", "status", "z_valve_start")
end


function set_start_all!(data::Dict{String,<:Any})
    apply_wm!(_set_start_all!, data)
end


function _set_start_all!(data::Dict{String,<:Any})
    _set_flow_start!(data)
    _set_head_start!(data)
    _set_indicator_start!(data)
    _set_flow_direction_start!(data)
end


function fix_all_flow_directions!(data::Dict{String,<:Any})
    apply_wm!(_fix_all_flow_directions!, data)
end


function _fix_all_flow_directions!(data::Dict{String,<:Any})
    _fix_flow_directions!(data, "pipe")
    _fix_flow_directions!(data, "short_pipe")
    _fix_flow_directions!(data, "pump")
    _fix_flow_directions!(data, "regulator")
    _fix_flow_directions!(data, "valve")
end


function _fix_flow_directions!(data::Dict{String,<:Any}, component_type::String)
    comps = values(data[component_type])
    _fix_flow_direction!.(comps)
end


function _fix_flow_direction!(component::Dict{String,<:Any})
    if haskey(component, "q") && !isapprox(component["q"], 0.0; atol=1.0e-6)
        component["y_min"] = component["q"] > 0.0 ? 1.0 : 0.0
        component["y_max"] = component["q"] > 0.0 ? 1.0 : 0.0
    end
end


function fix_all_indicators!(data::Dict{String,<:Any})
    apply_wm!(_fix_all_indicators!, data; apply_to_subnetworks = true)
end


function _fix_all_indicators!(data::Dict{String,<:Any})
    _fix_indicators!(data, "pump")
    _fix_indicators!(data, "regulator")
    _fix_indicators!(data, "valve")
end


function _fix_indicators!(data::Dict{String,<:Any}, component_type::String)
    comps = values(data[component_type])
    _fix_indicator!.(comps)
end


function _fix_indicator!(component::Dict{String,<:Any})
    if haskey(component, "status")
        _correct_status!(component)
        component["z_min"] = component["status"] === STATUS_ACTIVE ? 1.0 : 0.0
        component["z_max"] = component["status"] === STATUS_ACTIVE ? 1.0 : 0.0
    end
end


function turn_on_all_components!(data::Dict{String,<:Any})
    apply_wm!(_turn_on_all_components!, data; apply_to_subnetworks = true)
end


function _turn_on_all_components!(data::Dict{String,<:Any})
    _turn_on_components!(data, "pump")
    _turn_on_components!(data, "regulator")
    _turn_on_components!(data, "valve")
end


function _turn_on_components!(data::Dict{String,<:Any}, component_type::String)
    comps = values(data[component_type])
    _turn_on_component!.(comps)
end


function _turn_on_component!(component::Dict{String,<:Any})
    component["status"] = STATUS_ACTIVE
end


function relax_network!(data::Dict{String,<:Any})
    apply_wm!(_relax_network!, data; apply_to_subnetworks = true)
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
    apply_wm!(_recompute_bounds!, data)
end


function _recompute_bounds!(data::Dict{String, <:Any})
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
    haskey(comp, key) && (comp[key] = transform_flow.(comp[key]))
end


function _transform_flows!(comp::Dict{String,<:Any}, transform_flow::Function)
    for flow_key in ["flow_min", "flow_max", "flow_nominal", "flow_min_forward",
        "flow_max_reverse", "minor_loss", "flow_partition", "q", "qp", "qn"]
        _transform_flow_key!(comp, flow_key, transform_flow)
    end
end


function _transform_head_difference_key!(comp::Dict{String,<:Any}, key::String, transform_head::Function)
    haskey(comp, key) && (comp[key] = transform_head.(comp[key]))
end


function _transform_head_differences!(comp::Dict{String,<:Any}, transform_head::Function)
    for head_difference_key in ["dh", "dhp", "dhn", "g"]
        _transform_head_difference_key!(comp, head_difference_key, transform_head)
    end
end


function _apply_node_unit_transform!(data::Dict{String,<:Any}, transform_head::Function)
    wm_data = get_wm_data(data)

    for node in values(wm_data["node"])
        for head_key in ["elevation", "head", "head_min", "head_max", "head_nominal", "h", "p"]
            haskey(node, head_key) && (node[head_key] = transform_head.(node[head_key]))
        end
    end

    if haskey(wm_data, "time_series") && haskey(wm_data["time_series"], "node")
        for node in values(wm_data["time_series"]["node"])
            node["head_min"] = transform_head.(node["head_min"])
            node["head_max"] = transform_head.(node["head_max"])
            node["head_nominal"] = transform_head.(node["head_nominal"])
        end
    end
end


function _apply_reservoir_unit_transform!(data::Dict{String,<:Any}, transform_head::Function, transform_flow::Function)
    wm_data = get_wm_data(data)

    if !haskey(wm_data, "reservoir")
        return
    end

    for (_, reservoir) in wm_data["reservoir"]
        for head_key in ["head_min", "head_max", "head_nominal"]
            if haskey(reservoir, head_key)
                reservoir[head_key] = transform_head(reservoir[head_key])
            end
        end

        if haskey(reservoir, "q")
            reservoir["q"] = transform_flow(reservoir["q"])
        end
    end
end


function _apply_demand_unit_transform!(data::Dict{String,<:Any}, transform_flow::Function)
    wm_data = get_wm_data(data)

    if !haskey(wm_data, "demand")
        return
    end

    for (_, demand) in wm_data["demand"]
        for flow_key in ["flow_min", "flow_max", "flow_nominal", "q"]
            if haskey(demand, flow_key)
                demand[flow_key] = transform_flow(demand[flow_key])
            end
        end
    end

    if haskey(wm_data, "time_series") && haskey(wm_data["time_series"], "demand")
        for demand in values(wm_data["time_series"]["demand"])
            _transform_flows!(demand, transform_flow)
        end
    end
end


function _apply_tank_unit_transform!(data::Dict{String,<:Any}, transform_length::Function, transform_flow::Function)
    wm_data = get_wm_data(data)

    if !haskey(wm_data, "tank")
        return
    end

    for tank in values(wm_data["tank"])
        for length_key in ["min_level", "max_level", "init_level", "diameter"]
            if haskey(tank, length_key)
                tank[length_key] = transform_length(tank[length_key])
            end
        end

        for volume_key in ["min_vol", "V"]
            if haskey(tank, volume_key)
                tank[volume_key] *= transform_length(1.0)^3
            end
        end

        if haskey(tank, "q")
            tank["q"] = transform_flow(tank["q"])
        end
    end
end


function _make_per_unit_heads!(data::Dict{String,<:Any}, transform_head::Function)
    _apply_node_unit_transform!(get_wm_data(data), transform_head)
end


function _apply_link_unit_transform!(data::Dict{String,<:Any}, transform_flow::Function, transform_head::Function)
    wm_data = get_wm_data(data)

    for type in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
        if haskey(wm_data, type)
            map(x -> _transform_flows!(x, transform_flow), values(wm_data[type]))
            map(x -> _transform_head_differences!(x, transform_head), values(wm_data[type]))
        end

        if haskey(wm_data, "time_series") && haskey(wm_data["time_series"], type)
            _transform_flows!(wm_data["time_series"][type], transform_flow)
        end
    end
end


function _apply_pipe_unit_transform!(data::Dict{String,<:Any}, transform_length::Function, head_loss::String)
    wm_data = get_wm_data(data)
    
    if !haskey(wm_data, "pipe")
        return
    end

    for (_, pipe) in wm_data["pipe"]
        if haskey(pipe, "length")
            pipe["length"] = transform_length(pipe["length"])
        end

        if haskey(pipe, "diameter")
            pipe["diameter"] = transform_length(pipe["diameter"])
        end

        if uppercase(head_loss) == "D-W" && haskey(pipe, "roughness")
            pipe["roughness"] = transform_length(pipe["roughness"])
        end
    end
end


function _apply_des_pipe_unit_transform!(data::Dict{String,<:Any}, transform_length::Function, head_loss::String)
    wm_data = get_wm_data(data)

    if !haskey(wm_data, "des_pipe")
        return
    end

    for (_, des_pipe) in wm_data["des_pipe"]
        if haskey(des_pipe, "length")
            des_pipe["length"] = transform_length(des_pipe["length"])
        end

        if haskey(des_pipe, "diameter")
            des_pipe["diameter"] = transform_length(des_pipe["diameter"])
        end

        if uppercase(head_loss) == "D-W" && haskey(des_pipe, "roughness")
            des_pipe["roughness"] = transform_length(des_pipe["roughness"])
        end
    end
end


function _apply_pump_unit_transform!(
    data::Dict{String,<:Any}, transform_mass::Function, transform_flow::Function,
    transform_length::Function, transform_time::Function)
    wm_data = get_wm_data(data)

    if !haskey(wm_data, "pump")
        return
    end

    power_scalar = transform_mass(1.0) * transform_length(1.0)^2 / transform_time(1.0)^3
    energy_scalar = transform_mass(1.0) * transform_length(1.0)^2 / transform_time(1.0)^2

    for (_, pump) in wm_data["pump"]
        if haskey(pump, "head_curve")
            pump["head_curve"] = [(transform_flow(x[1]), x[2]) for x in pump["head_curve"]]
            pump["head_curve"] = [(x[1], transform_length(x[2])) for x in pump["head_curve"]]
        end

        if haskey(pump, "efficiency_curve")
            pump["efficiency_curve"] = [(transform_flow(x[1]), x[2]) for x in pump["efficiency_curve"]]
        end

        if haskey(pump, "energy_price")
            pump["energy_price"] /= energy_scalar
        end

        if haskey(pump, "E")
            pump["E"] *= energy_scalar
        end

        if haskey(pump, "power_fixed")
            pump["power_fixed"] *= power_scalar
        end

        if haskey(pump, "P")
            pump["P"] *= power_scalar
        end

        if haskey(pump, "c")
            pump["c"] *= energy_scalar
        end

        if haskey(pump, "min_inactive_time")
            pump["min_inactive_time"] = transform_time(pump["min_inactive_time"])
        end 

        if haskey(pump, "min_active_time")
            pump["min_active_time"] = transform_time(pump["min_active_time"])
        end 

        if haskey(pump, "power_per_unit_flow")
            pump["power_per_unit_flow"] *= power_scalar / transform_flow(1.0)
        end
    end

    if haskey(wm_data, "time_series") && haskey(wm_data["time_series"], "pump")
        for pump in values(wm_data["time_series"]["pump"])
            if haskey(pump, "energy_price")
                pump["energy_price"] ./= energy_scalar
            end

            if haskey(pump, "power_fixed")
                pump["power_fixed"] .*= power_scalar
            end

            if haskey(pump, "power_per_unit_flow")
                pump["power_per_unit_flow"] .*= power_scalar / transform_flow(1.0)
            end
        end
    end
end


function _apply_regulator_unit_transform!(data::Dict{String,<:Any}, transform_head::Function, transform_length::Function)
    wm_data = get_wm_data(data)

    if !haskey(wm_data, "regulator")
        return
    end

    for (i, regulator) in wm_data["regulator"]
        if haskey(regulator, "diameter")
            regulator["diameter"] = transform_length(regulator["diameter"])
        end

        if haskey(regulator, "setting")
            regulator["setting"] = transform_head(regulator["setting"])
        end
    end
end


function _calc_scaled_gravity(base_length::Float64, base_time::Float64)
    return _GRAVITY * (base_time)^2 / base_length
end


function _calc_scaled_gravity(data::Dict{String, <:Any})
    wm_data = get_wm_data(data)

    if wm_data["per_unit"]
        base_time = 1.0 / _calc_time_per_unit_transform(wm_data)(1.0)
        base_length = 1.0 / _calc_length_per_unit_transform(wm_data)(1.0) 
        return _calc_scaled_gravity(base_length, base_time)
    else
        return _GRAVITY
    end
end


function _calc_scaled_density(base_mass::Float64, base_length::Float64)
    return _DENSITY * base_length^3 / base_mass
end


function _calc_scaled_density(data::Dict{String, <:Any})
    wm_data = get_wm_data(data)

    if wm_data["per_unit"]
        base_mass = 1.0 / _calc_mass_per_unit_transform(wm_data)(1.0)
        base_length = 1.0 / _calc_length_per_unit_transform(wm_data)(1.0)
        return _calc_scaled_density(base_mass, base_length)
    else
        return _DENSITY
    end
end

"Transforms data to a per-unit system."
function make_per_unit!(data::Dict{String,<:Any})
    apply_wm!(_make_per_unit!, data; apply_to_subnetworks = false)
end


function _make_per_unit!(data::Dict{String,<:Any})
    if get(data, "per_unit", false) == false
        _si_to_per_unit_metadata!(data)

        if _IM.ismultinetwork(data)
            for (i, _) in data["nw"]
                _si_to_per_unit_nw!(data, id = i)
            end
        else
            _si_to_per_unit_nw!(data)
        end

        wm_data = get_wm_data(data)
        wm_data["per_unit"] = true
    end
end


function _si_to_per_unit_nw!(data::Dict{String,<:Any}; id = "0")
    # Precompute SI-to-per-unit transformation functions.
    mass_transform = _calc_mass_per_unit_transform(data)
    time_transform = _calc_time_per_unit_transform(data)
    length_transform = _calc_length_per_unit_transform(data)
    head_transform = _calc_head_per_unit_transform(data)
    flow_transform = _calc_flow_per_unit_transform(data)

    # Perform the transformation.
    _transform_nw!(data, mass_transform, time_transform, length_transform,
        head_transform, flow_transform; to_per_unit = true, id = id)
end


function _si_to_per_unit_metadata!(data::Dict{String,<:Any})
    # Precompute SI-to-per-unit transformation functions.
    mass_transform = _calc_mass_per_unit_transform(data)
    time_transform = _calc_time_per_unit_transform(data)
    length_transform = _calc_length_per_unit_transform(data)
    head_transform = _calc_head_per_unit_transform(data)
    flow_transform = _calc_flow_per_unit_transform(data)

    # Perform the transformation.
    _transform_metadata!(data, mass_transform, time_transform, length_transform,
        head_transform, flow_transform, true)
end


"Transforms data to SI units."
function make_si_units!(data::Dict{String,<:Any})
    apply_wm!(_make_si_units!, data; apply_to_subnetworks = false)
end


function _make_si_units!(data::Dict{String,<:Any})
    if get(data, "per_unit", false) == true
        _per_unit_to_si_metadata!(data)

        if _IM.ismultinetwork(data)
            for (i, _) in data["nw"]
                _per_unit_to_si_nw!(data, id = i)
            end
        else
            _per_unit_to_si_nw!(data)
        end

        wm_data = get_wm_data(data)
        wm_data["per_unit"] = false
    end
end


function _per_unit_to_si_nw!(data::Dict{String,<:Any}; id = "0")
    # Precompute per-unit-to-SI transformation functions.
    mass_transform = _calc_mass_si_transform(data)
    time_transform = _calc_time_si_transform(data)
    length_transform = _calc_length_si_transform(data)
    head_transform = _calc_head_si_transform(data)
    flow_transform = _calc_flow_si_transform(data)

    # Perform the transformation.
    _transform_nw!(data, mass_transform, time_transform, length_transform,
        head_transform, flow_transform; to_per_unit = false, id = id)
end


function _per_unit_to_si_metadata!(data::Dict{String,<:Any})
    # Precompute per-unit-to-SI transformation functions.
    mass_transform = _calc_mass_si_transform(data)
    time_transform = _calc_time_si_transform(data)
    length_transform = _calc_length_si_transform(data)
    head_transform = _calc_head_si_transform(data)
    flow_transform = _calc_flow_si_transform(data)

    # Perform the transformation.
    _transform_metadata!(data, mass_transform, time_transform, length_transform,
        head_transform, flow_transform, false)
end


function _transform_metadata!(
    data::Dict{String,<:Any},
    mass_transform::Function,
    time_transform::Function,
    length_transform::Function,
    head_transform::Function,
    flow_transform::Function,
    to_per_unit::Bool,
)
    # Get the WaterModels subdictionary.
    wm_data = get_wm_data(data)

    if haskey(wm_data, "viscosity")
        # Convert viscosity to per-unit units.
        wm_data["viscosity"] *= mass_transform(1.0) /
            (length_transform(1.0) * time_transform(1.0))
    end

    # Store important base values for later conversions.
    if to_per_unit && (!haskey(wm_data, "per_unit") || !wm_data["per_unit"])
        wm_data["base_head"] = get(wm_data, "base_head", 1.0 / head_transform(1.0))
        wm_data["base_mass"] = get(wm_data, "base_mass", 1.0 / mass_transform(1.0))
        wm_data["base_length"] = get(wm_data, "base_length", 1.0 / length_transform(1.0))
        wm_data["base_time"] = get(wm_data, "base_time", 1.0 / time_transform(1.0))
        wm_data["base_flow"] = get(wm_data, "base_flow", 1.0 / flow_transform(1.0))
    end

    # Convert time-related data to SI units.
    if haskey(wm_data, "time_series")
        time_step = time_transform(wm_data["time_series"]["time_step"])
        wm_data["time_series"]["time_step"] = time_step
        duration = time_transform(data["time_series"]["duration"])
        wm_data["time_series"]["duration"] = duration
    end

    if haskey(wm_data, "time_step")
        wm_data["time_step"] = time_transform(wm_data["time_step"])
    end
end


function _transform_nw!(
    data::Dict{String,<:Any},
    mass_transform::Function,
    time_transform::Function,
    length_transform::Function,
    head_transform::Function,
    flow_transform::Function;
    to_per_unit::Bool = false,
    id = "0"
)
    # Get the WaterModels subdictionary.
    wm_data = get_wm_data(data)

    # Get the head loss form.
    head_loss = wm_data["head_loss"]
    @assert uppercase(head_loss) in ["D-W", "H-W"]

    # Ensure the transformation is taking place in the correct direction.
    @assert to_per_unit != get(wm_data, "per_unit", false)

    # Get the corresponding subnetwork data, if necessary.
    wm_nw_data = (id == "0") ? wm_data : wm_data["nw"][id]

    # Apply transformations to node-connecting components.
    _apply_link_unit_transform!(wm_nw_data, flow_transform, head_transform)
    _apply_pipe_unit_transform!(wm_nw_data, length_transform, head_loss)
    _apply_des_pipe_unit_transform!(wm_nw_data, length_transform, head_loss)
    _apply_pump_unit_transform!(wm_nw_data, mass_transform,
        flow_transform, length_transform, time_transform)
    _apply_regulator_unit_transform!(wm_nw_data, head_transform, length_transform)

    # Apply transformations to nodal components.
    _apply_node_unit_transform!(wm_nw_data, head_transform)
    _apply_demand_unit_transform!(wm_nw_data, flow_transform)
    _apply_tank_unit_transform!(wm_nw_data, length_transform, flow_transform)
    _apply_reservoir_unit_transform!(wm_nw_data, head_transform, flow_transform)
end


function set_warm_start!(data::Dict{String,<:Any})
    apply_wm!(_set_warm_start!, data)
end


function _set_warm_start!(data::Dict{String, <:Any})
    _set_node_warm_start!(data)
    _set_demand_warm_start!(data)
    _set_reservoir_warm_start!(data)
    _set_tank_warm_start!(data)

    _set_pipe_warm_start!(data)
    _set_pump_warm_start!(data)
    _set_short_pipe_warm_start!(data)
    _set_valve_warm_start!(data)
end