# Functions for working with WaterModels data elements.

"WaterModels wrapper for the InfrastructureModels `apply!` function."
function apply_wm!(func!::Function, data::Dict{String, <:Any}; apply_to_subnetworks::Bool = true)
    _IM.apply!(func!, data, wm_it_name; apply_to_subnetworks = apply_to_subnetworks)
end


function update_component_keys!(data_1::Dict{String,<:Any}, data_2::Dict{String,<:Any})
    for (n, data_nw) in data_1["nw"]
        if haskey(data_2["nw"], n)
             component_types = keys(data_2["nw"][n])

            for component_type in component_types
                mod_comp = data_2["nw"][n][component_type]
                name_map = Dict(x["source_id"][2] => string(x["index"])
                    for (i, x) in data_nw[component_type])

                for (name, index) in name_map
                    if haskey(mod_comp, name)
                        mod_comp[index] = pop!(mod_comp, name)
                    end
                end
            end
        end
    end
end


function update_data_from_source!(data::Dict{String,<:Any}, modifications::Dict{String,<:Any})
    if ismultinetwork(data) && ismultinetwork(modifications)
        for (n, data_nw) in data["nw"]
            if haskey(modifications["nw"], n)
                component_types = keys(modifications["nw"][n])

                for component_type in component_types
                    mod_comp = modifications["nw"][n][component_type]
                    name_map = Dict(x["source_id"][2] => string(x["index"])
                        for (i, x) in data_nw[component_type])

                    for (name, index) in name_map
                        mod_comp[index] = pop!(mod_comp, name)
                    end
                end
            end
        end
    end

    _IM.update_data!(data, modifications)
end


"Check that all nodes are unique and other components link to valid nodes."
function check_connectivity(data::Dict{String,<:Any})
    apply_wm!(_check_connectivity, data)
end


"Check that all nodes are unique and other components link to valid nodes for a
single-network data set that does not contain other multi-infrastructure data."
function _check_connectivity(data::Dict{String,<:Any})
    node_ids = Set(node["index"] for (i, node) in data["node"])
    @assert(length(node_ids) == length(data["node"]))

    for comp_type in _NODE_CONNECTED_COMPONENTS
        for (i, comp) in data[comp_type]
            if !(comp["node"] in node_ids)
                error_message = "Node $(comp["node"]) in $(comp_type) $(i) is not defined."
                Memento.error(_LOGGER, error_message)
            end
        end
    end

    for comp_type in _LINK_COMPONENTS
        for (i, comp) in data[comp_type]
            if !(comp["node_fr"] in node_ids)
                error_message = "From node $(comp["node_fr"]) in "
                error_message *= "$(replace(comp_type, "_" => " ")) $(i) is not defined."
                Memento.error(_LOGGER, error_message)
            end

            if !(comp["node_to"] in node_ids)
                error_message = "To node $(comp["node_to"]) in "
                error_message *= "$(replace(comp_type, "_" => " ")) $(i) is not defined."
                Memento.error(_LOGGER, error_message)
            end
        end
    end
end


"Check that active components are not connected to inactive nodes."
function check_status(data::Dict{String,<:Any})
    apply_wm!(_check_status, data)
end


"Check that active components are not connected to inactive nodes for a single-network data
set that does not contain other multi-infrastructure data."
function _check_status(data::Dict{String,<:Any})
    active_nodes = filter(x -> x.second["status"] != STATUS_INACTIVE, data["node"])
    active_node_ids = Set(node["index"] for (i, node) in active_nodes)

    for comp_type in _NODE_CONNECTED_COMPONENTS
        for (i, comp) in data[comp_type]
            if comp["status"] != STATUS_INACTIVE && !(comp["node"] in active_node_ids)
                warning_message = "Active $(comp_type) $(i) is connected to inactive "
                warning_message *= "node $(comp["node"])."
                Memento.warn(_LOGGER, warning_message)
            end
        end
    end

    for comp_type in _LINK_COMPONENTS
        for (i, comp) in data[comp_type]
            if comp["status"] != STATUS_INACTIVE && !(comp["node_fr"] in active_node_ids)
                warning_message = "Active $(comp_type) $(i) is connected to inactive "
                warning_message *= "from node $(comp["node_fr"])."
                Memento.warn(_LOGGER, warning_message)
            end

            if comp["status"] != STATUS_INACTIVE && !(comp["node_to"] in active_node_ids)
                warning_message = "Active $(comp_type) $(i) is connected to inactive "
                warning_message *= "to node $(comp["node_to"])."
                Memento.warn(_LOGGER, warning_message)
            end
        end
    end
end


function correct_enums!(data::Dict{String,<:Any})
    correct_statuses!(data)
    correct_flow_directions!(data)
    correct_pump_types!(data)
end


function correct_flow_directions!(data::Dict{String,<:Any})
    apply_wm!(_correct_flow_directions!, data; apply_to_subnetworks = true)
end


function _correct_flow_directions!(data::Dict{String,<:Any})
    for component_type in _LINK_COMPONENTS
        components = values(get(data, component_type, Dict{String,Any}()))
        _correct_flow_direction!.(components)
    end
end


"Correct flow direction attribute of link-type components."
function _correct_flow_direction!(comp::Dict{String, <:Any})
    flow_direction = get(comp, "flow_direction", FLOW_DIRECTION_UNKNOWN)

    if isa(flow_direction, FLOW_DIRECTION)
        # Use the flow direction from above, which is correctly typed.
        comp["flow_direction"] = flow_direction
    else
        # Correct "flow_direction" type to the enum type.
        comp["flow_direction"] = FLOW_DIRECTION(flow_direction)
    end

    if get(comp, "flow_min", -Inf) > 0.0
        # If minimum flow is positive, assume positively-directed flow.
        comp["flow_direction"] = FLOW_DIRECTION_POSITIVE
    elseif get(comp, "flow_max", Inf) < 0.0
        # If maximum flow is negative, assume negatively-directed flow.
        comp["flow_direction"] = FLOW_DIRECTION_NEGATIVE
    end
end


function correct_statuses!(data::Dict{String,<:Any})
    apply_wm!(_correct_statuses!, data; apply_to_subnetworks = true)
end


function _correct_statuses!(data::Dict{String,<:Any})
    for component_type in vcat(_LINK_COMPONENTS, _NODE_CONNECTED_COMPONENTS)
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

    for pump in values(get(data, "ne_pump", Dict{String, Any}()))
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
    base_flow = get(wm_data, "base_flow", 1.0)
    base_length = get(wm_data, "base_length", 1.0)
    base_mass = get(wm_data, "base_mass", 1.0)
    base_time = get(wm_data, "base_time", 1.0)

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

    # Set partitions for all expansion pumps in the network.
    for pump in values(get(data, "ne_pump", Dict{String, Any}()))
        set_pump_flow_partition!(pump, error_tolerance, length_tolerance)
    end
end


"Transform length values in SI units to per-unit units."
function _calc_length_per_unit_transform(data::Dict{String,<:Any})
    wm_data = get_wm_data(data)

    if haskey(wm_data, "base_length")
        return x -> x / wm_data["base_length"]
    else
        median_midpoint = _calc_node_head_median_midpoint(wm_data)
        return x -> x / median_midpoint
    end
end


"Correct status attribute of a component."
function _correct_status!(comp::Dict{String, <:Any})
    status = get(comp, "status", STATUS_UNKNOWN)

    if isa(status, STATUS)
        comp["status"] = status
    else
        comp["status"] = STATUS(Int(round(status)))
    end
end


"Transform length values SI units to per-unit units."
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


"Transform head values in per-unit units to SI units."
function _calc_head_per_unit_untransform(data::Dict{String,<:Any})
    median_midpoint = _calc_node_head_median_midpoint(get_wm_data(data))
    return x -> x * median_midpoint + median_midpoint
end


"Transform time values in SI units to per-unit units."
function _calc_time_per_unit_transform(data::Dict{String,<:Any})
    wm_data = get_wm_data(data)

    if haskey(wm_data, "base_time")
        @assert wm_data["base_time"] > 0.0
        return x -> x / wm_data["base_time"]
    else
         # Translation: number of meters per WaterModels-length
        head_midpoint = _calc_node_head_median_midpoint(wm_data)
        flow_midpoint = _calc_median_abs_flow_midpoint(wm_data)

        # Translation: convert number of seconds to WaterModels-time.
        return x -> x / (head_midpoint^3 / flow_midpoint)
    end
end

"Transform flow values in SI units to per-unit units."
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


"Transform mass values in SI units to per-unit units."
function _calc_mass_per_unit_transform(data::Dict{String,<:Any})
    # Translation: convert kilograms to WaterModels-mass.
    return x -> x / get(get_wm_data(data), "base_mass", 1.0e19)
end


function _calc_median_abs_flow_midpoint(data::Dict{String,<:Any})
    wm_data = get_wm_data(data)

    # Compute flow midpoints for all possible node-connecting components.
    q_des_pipe = [_calc_abs_flow_midpoint(x) for (i, x) in wm_data["des_pipe"]]
    q_pipe = [_calc_abs_flow_midpoint(x) for (i, x) in wm_data["pipe"]]
    q_pump = [_calc_abs_flow_midpoint(x) for (i, x) in wm_data["pump"]]
    q_ne_pump = [_calc_abs_flow_midpoint(x) for (i, x) in wm_data["ne_pump"]]
    q_regulator = [_calc_abs_flow_midpoint(x) for (i, x) in wm_data["regulator"]]
    q_short_pipe = [_calc_abs_flow_midpoint(x) for (i, x) in wm_data["short_pipe"]]
    q_ne_short_pipe = [_calc_abs_flow_midpoint(x) for (i, x) in wm_data["ne_short_pipe"]]
    q_valve = [_calc_abs_flow_midpoint(x) for (i, x) in wm_data["valve"]]
    q = vcat(q_des_pipe, q_pipe, q_pump, q_ne_pump, q_regulator, q_short_pipe, q_ne_short_pipe, q_valve)

    # Calculate the median of all values computed above.
    q_median = Statistics.median(q)

    # Return the median if not equal to zero. Otherwise, return 1.0.
    return q_median != 0.0 ? q_median : 1.0
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


"Compute the maximum capacity of the network."
function calc_capacity_max(data::Dict{String, <:Any})
    # Get the WaterModels portion of the data dictionary.
    wm_data = get_wm_data(data)

    # Compute the sum of all maximal flows from demands.
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


function make_ts_metadata!(data::Dict{String, <:Any})
    @assert ismultinetwork(data) # Ensure data is multinetwork.

    if !haskey(data, "time_series")
        data["time_series"] = Dict{String, Any}()
    end

    data["time_series"]["duration"] = data["duration"]
    data["time_series"]["num_steps"] = length(data["nw"])
    data["time_series"]["time_step"] = data["time_step"]
end


function make_component_ts!(data::Dict{String, <:Any}, comp_type::String, key::String)
    @assert ismultinetwork(data) # Ensure data is multinetwork.
    nws = sort([parse(Int, x) for x in keys(data["nw"])])

    if !haskey(data, "time_series")
        data["time_series"] = Dict{String, Any}()
    end

    if !haskey(data["time_series"], comp_type)
        data["time_series"][comp_type] = Dict{String, Any}()
    end

    if !haskey(data["time_series"][comp_type], key)
        data["time_series"][comp_type][key] = Dict{String, Any}()
    end

    for val_key in keys(data["nw"][string(nws[1])][comp_type][key])
        vals = Array{Any}([])

        for nw in nws
            if haskey(data["nw"][string(nw)][comp_type][key], val_key)
                push!(vals, data["nw"][string(nw)][comp_type][key][val_key])
            else
                # Assume the last value.
                push!(vals, vals[end])
            end
        end

        data["time_series"][comp_type][key][val_key] = vals
     end
end


function make_single_network(data::Dict{String, <:Any})
    @assert ismultinetwork(data) # Ensure data is multinetwork.

    data_s = deepcopy(data)
    make_ts_metadata!(data_s)
    nws = sort([parse(Int, x) for x in keys(data["nw"])])
    nw_1_str = string(nws[1])

    for comp_type in ["tank", "regulator", "pump", "ne_pump", "des_pipe", "pump_group", "demand",
        "tank_group", "reservoir", "node", "short_pipe", "ne_short_pipe", "valve", "pipe"]
        if !haskey(data_s["nw"][nw_1_str], comp_type)
            continue
        end

        comp_keys = keys(data_s["nw"][nw_1_str][comp_type])
        make_component_ts!.(Ref(data_s), comp_type, comp_keys)

        if !haskey(data_s, comp_type)
            data_s[comp_type] = Dict{String, Any}()
        end

        for comp_key in comp_keys
            data_s[comp_type][comp_key] = pop!(data_s["nw"][nw_1_str][comp_type], comp_key)
        end
    end

    data_s["multinetwork"] = false
    delete!(data_s, "nw")

    return data_s
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
    set_start!(data, "ne_pump", "q", "q_ne_pump_start")
    set_start!(data, "regulator", "q", "q_regulator_start")
    set_start!(data, "short_pipe", "q", "q_short_pipe_start")
    set_start!(data, "ne_short_pipe", "q", "q_ne_short_pipe_start")
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
    set_direction_start_from_flow!(data, "ne_pump", "q", "y_ne_pump_start")
    set_direction_start_from_flow!(data, "regulator", "q", "y_regulator_start")
    set_direction_start_from_flow!(data, "short_pipe", "q", "y_short_pipe_start")
    set_direction_start_from_flow!(data, "ne_short_pipe", "q", "y_ne_short_pipe_start")
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
    _fix_flow_directions!(data, "ne_short_pipe")
    _fix_flow_directions!(data, "pump")
    _fix_flow_directions!(data, "ne_pump")
    _fix_flow_directions!(data, "regulator")
    _fix_flow_directions!(data, "valve")
end


function _fix_flow_directions!(data::Dict{String,<:Any}, component_type::String)
    comps = values(data[component_type])
    _fix_flow_direction!.(comps)
end


function _fix_flow_direction!(component::Dict{String,<:Any})
    if haskey(component, "q") && !isapprox(component["q"], 0.0; atol = 1.0e-6)
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


function set_bounds_from_time_series!(data::Dict{String,<:Any})
    apply_wm!(_set_bounds_from_time_series!, data; apply_to_subnetworks = false)
end


function _set_bounds_from_time_series!(data::Dict{String,<:Any})
    # Only operate on non-multinetworks.
    @assert !ismultinetwork(data)

    if haskey(data, "time_series")
        time_series = Ref(data["time_series"])

        # Set node and nodal component data based on time series bounds.
        _set_node_bounds_from_time_series!.(values(data["node"]), time_series)
        _set_demand_bounds_from_time_series!.(values(data["demand"]), time_series)
        _set_reservoir_bounds_from_time_series!.(values(data["reservoir"]), time_series)
        _set_tank_bounds_from_time_series!.(values(data["tank"]), time_series)

        for link in _LINK_COMPONENTS
            # Set link component data based on time series bounds.
            _set_link_bounds_from_time_series!.(values(data[link]), link, time_series)
        end
    end
end


function make_all_nondispatchable!(data::Dict{String,<:Any})
    apply_wm!(_make_all_nondispatchable!, data)
end


function _make_all_nondispatchable!(data::Dict{String,<:Any})
    for comp_type in vcat(_LINK_COMPONENTS, _NODE_CONNECTED_COMPONENTS)
        for comp in values(data[comp_type])
            if haskey(comp, "dispatchable")
                comp["dispatchable"] = false
            end
        end
    end
end


"Convenience function for recomputing component bounds, e.g., after modifying data."
function recompute_bounds!(data::Dict{String, <:Any})
    apply_wm!(_recompute_bounds!, data)
end


function _recompute_bounds!(data::Dict{String, <:Any})
    # Clear the existing flow bounds for node-connecting components.
    for comp_type in _LINK_COMPONENTS
        map(x -> x["flow_min"] = -Inf, values(data[comp_type]))
        map(x -> x["flow_max"] = Inf, values(data[comp_type]))
    end

    # Recompute bounds and correct data.
    correct_network_data!(data)
end


function _transform_flow_key!(comp::Dict{String,<:Any}, key::String, transform_flow::Function)
    haskey(comp, key) && (comp[key] = transform_flow.(comp[key]))
end


function _transform_flows!(comp::Dict{String,<:Any}, transform_flow::Function)
    _transform_flow_key!(comp, "flow_min", transform_flow)
    _transform_flow_key!(comp, "flow_min_forward", transform_flow)
    _transform_flow_key!(comp, "flow_max", transform_flow)
    _transform_flow_key!(comp, "flow_max_reverse", transform_flow)
    _transform_flow_key!(comp, "flow_nominal", transform_flow)
    _transform_flow_key!(comp, "minor_loss", transform_flow)
    _transform_flow_key!(comp, "flow_partition", transform_flow)
end


function _make_per_unit_nodes!(data::Dict{String,<:Any}, transform_head::Function)
    wm_data = get_wm_data(data)

    for node in values(wm_data["node"])
        node["elevation"] = transform_head(node["elevation"])
        node["head_min"] = transform_head(node["head_min"])
        node["head_max"] = transform_head(node["head_max"])
        node["head_nominal"] = transform_head(node["head_nominal"])
    end

    if haskey(wm_data, "time_series") && haskey(wm_data["time_series"], "node")
        for node in values(wm_data["time_series"]["node"])
            node["head_min"] = transform_head.(node["head_min"])
            node["head_max"] = transform_head.(node["head_max"])
            node["head_nominal"] = transform_head.(node["head_nominal"])
        end
    end
end


function _make_per_unit_reservoir!(data::Dict{String,<:Any}, transform_head::Function)
    wm_data = get_wm_data(data)

    for (i, reservoir) in wm_data["reservoir"]
        if haskey(reservoir, "head_min")
            reservoir["head_min"] = transform_head(reservoir["head_min"])
        end

        if haskey(reservoir, "head_max")
            reservoir["head_max"] = transform_head(reservoir["head_max"])
        end

        if haskey(reservoir, "head_nominal")
            reservoir["head_nominal"] = transform_head(reservoir["head_nominal"])
        end
    end
end


function _make_per_unit_demands!(data::Dict{String,<:Any}, transform_flow::Function)
    wm_data = get_wm_data(data)

    for (i, demand) in wm_data["demand"]
        demand["flow_min"] = transform_flow(demand["flow_min"])
        demand["flow_max"] = transform_flow(demand["flow_max"])
        demand["flow_nominal"] = transform_flow(demand["flow_nominal"])
    end

    if haskey(wm_data, "time_series") && haskey(wm_data["time_series"], "demand")
        for demand in values(wm_data["time_series"]["demand"])
            _transform_flows!(demand, transform_flow)
        end
    end
end


function _make_per_unit_tanks!(data::Dict{String,<:Any}, transform_length::Function)
    wm_data = get_wm_data(data)

    for tank in values(wm_data["tank"])
        tank["min_level"] = transform_length(tank["min_level"])
        tank["max_level"] = transform_length(tank["max_level"])
        tank["init_level"] = transform_length(tank["init_level"])
        tank["diameter"] = transform_length(tank["diameter"])
        tank["min_vol"] *= transform_length(1.0)^3
    end
end


function _make_per_unit_heads!(data::Dict{String,<:Any}, transform_head::Function)
    _make_per_unit_nodes!(get_wm_data(data), transform_head)
end


function _make_per_unit_flows!(data::Dict{String,<:Any}, transform_flow::Function)
    wm_data = get_wm_data(data)

    for type in ["des_pipe", "pipe", "pump", "ne_pump", "regulator", "short_pipe", "ne_short_pipe", "valve"]
        map(x -> _transform_flows!(x, transform_flow), values(wm_data[type]))

        if haskey(wm_data, "time_series") && haskey(wm_data["time_series"], type)
            _transform_flows!(wm_data["time_series"][type], transform_flow)
        end
    end
end


function _make_per_unit_pipes!(data::Dict{String,<:Any}, transform_length::Function)
    wm_data = get_wm_data(data)

    for (i, pipe) in wm_data["pipe"]
        pipe["length"] = transform_length(pipe["length"])
        pipe["diameter"] = transform_length(pipe["diameter"])

        if uppercase(wm_data["head_loss"]) == "D-W"
            pipe["roughness"] = transform_length(pipe["roughness"])
        end
    end
end


function _make_per_unit_des_pipes!(data::Dict{String,<:Any}, transform_length::Function)
    wm_data = get_wm_data(data)

    for (i, des_pipe) in wm_data["des_pipe"]
        des_pipe["length"] = transform_length(des_pipe["length"])
        des_pipe["diameter"] = transform_length(des_pipe["diameter"])
    end

    if uppercase(wm_data["head_loss"]) == "D-W"
        for (i, des_pipe) in wm_data["des_pipe"]
            des_pipe["roughness"] = transform_length(des_pipe["roughness"])
        end
    end
end


function _make_per_unit_pumps!(
    data::Dict{String,<:Any}, transform_mass::Function, transform_flow::Function,
    transform_length::Function, transform_time::Function)
    wm_data = get_wm_data(data)

    power_scalar = transform_mass(1.0) * transform_length(1.0)^2 / transform_time(1.0)^3
    energy_scalar = transform_mass(1.0) * transform_length(1.0)^2 / transform_time(1.0)^2

    for (i, pump) in wm_data["pump"]
        pump["head_curve"] = [(transform_flow(x[1]), x[2]) for x in pump["head_curve"]]
        pump["head_curve"] = [(x[1], transform_length(x[2])) for x in pump["head_curve"]]

        if haskey(pump, "efficiency_curve")
            pump["efficiency_curve"] = [(transform_flow(x[1]), x[2]) for x in pump["efficiency_curve"]]
        end

        if haskey(pump, "energy_price")
            pump["energy_price"] /= energy_scalar
        end

        if haskey(pump, "power_fixed")
            pump["power_fixed"] *= power_scalar
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

function _make_per_unit_ne_pumps!(
    data::Dict{String,<:Any}, transform_mass::Function, transform_flow::Function,
    transform_length::Function, transform_time::Function)
    wm_data = get_wm_data(data)

    power_scalar = transform_mass(1.0) * transform_length(1.0)^2 / transform_time(1.0)^3
    energy_scalar = transform_mass(1.0) * transform_length(1.0)^2 / transform_time(1.0)^2

    for (i, pump) in wm_data["ne_pump"]
        pump["head_curve"] = [(transform_flow(x[1]), x[2]) for x in pump["head_curve"]]
        pump["head_curve"] = [(x[1], transform_length(x[2])) for x in pump["head_curve"]]

        if haskey(pump, "efficiency_curve")
            pump["efficiency_curve"] = [(transform_flow(x[1]), x[2]) for x in pump["efficiency_curve"]]
        end

        if haskey(pump, "energy_price")
            pump["energy_price"] /= energy_scalar
        end

        if haskey(pump, "power_fixed")
            pump["power_fixed"] *= power_scalar
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

    if haskey(wm_data, "time_series") && haskey(wm_data["time_series"], "ne_pump")
        for pump in values(wm_data["time_series"]["ne_pump"])
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


function _make_per_unit_regulators!(data::Dict{String,<:Any}, transform_length::Function)
    wm_data = get_wm_data(data)

    for (i, regulator) in wm_data["regulator"]
        regulator["setting"] = transform_length(regulator["setting"])
    end
end


function _calc_scaled_gravity(base_length::Float64, base_time::Float64)
    return _GRAVITY * (base_time)^2 / base_length
end


function _calc_scaled_gravity(data::Dict{String, <:Any})
    wm_data = get_wm_data(data)
    base_time = 1.0 / _calc_time_per_unit_transform(wm_data)(1.0)
    base_length = 1.0 / _calc_length_per_unit_transform(wm_data)(1.0)
    return _GRAVITY * (base_time)^2 / base_length
end


function _calc_scaled_density(base_mass::Float64, base_length::Float64)
    return _DENSITY * base_length^3 / base_mass
end


function _calc_scaled_density(data::Dict{String, <:Any})
    wm_data = get_wm_data(data)
    base_mass = 1.0 / _calc_mass_per_unit_transform(wm_data)(1.0)
    base_length = 1.0 / _calc_length_per_unit_transform(wm_data)(1.0)
    return _DENSITY * base_length^3 / base_mass
end


function make_per_unit!(data::Dict{String,<:Any})
    apply_wm!(_make_per_unit!, data)
end


function _make_per_unit!(data::Dict{String,<:Any})
    if get(data, "per_unit", false) == false
        # Precompute per-unit transformation functions.
        mass_transform = _calc_mass_per_unit_transform(data)
        time_transform = _calc_time_per_unit_transform(data)
        length_transform = _calc_length_per_unit_transform(data)
        head_transform = _calc_head_per_unit_transform(data)
        flow_transform = _calc_flow_per_unit_transform(data)

        # Apply per-unit transformations to node-connecting components.
        _make_per_unit_flows!(data, flow_transform)
        _make_per_unit_pipes!(data, length_transform)
        _make_per_unit_des_pipes!(data, length_transform)
        _make_per_unit_pumps!(data, mass_transform,
            flow_transform, length_transform, time_transform)
        _make_per_unit_ne_pumps!(data, mass_transform,
            flow_transform, length_transform, time_transform)
        _make_per_unit_regulators!(data, length_transform)

        # Apply per-unit transformations to nodal components.
        _make_per_unit_nodes!(data, head_transform)
        _make_per_unit_demands!(data, flow_transform)
        _make_per_unit_tanks!(data, length_transform)
        _make_per_unit_reservoir!(data, head_transform)

        # Convert viscosity to per-unit units.
        data["viscosity"] *= mass_transform(1.0) /
            (length_transform(1.0) * time_transform(1.0))

        # Store important base values for later conversions.
        data["base_head"] = 1.0 / head_transform(1.0)
        data["base_mass"] = 1.0 / mass_transform(1.0)
        data["base_length"] = 1.0 / length_transform(1.0)
        data["base_time"] = 1.0 / time_transform(1.0)
        data["base_flow"] = 1.0 / flow_transform(1.0)

        if haskey(data, "time_series")
            time_step = time_transform(data["time_series"]["time_step"])
            data["time_series"]["time_step"] = time_step
            duration = time_transform(data["time_series"]["duration"])
            data["time_series"]["duration"] = duration
        end

        # Set the per-unit flag.
        data["time_step"] = time_transform(data["time_step"])
        data["per_unit"] = true
    end
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
    _set_ne_pump_warm_start!(data)
    _set_short_pipe_warm_start!(data)
    _set_ne_short_pipe_warm_start!(data)
    _set_valve_warm_start!(data)
end


"""
Deactivates components that are not needed in the network by repeated calls to
`propagate_topology_status!`. This implementation has quadratic complexity.
"""
function simplify_network!(data::Dict{String,<:Any})::Bool
    revised, num_iterations = true, 0

    while revised
        revised = false
        revised |= propagate_topology_status!(data)
        # revised |= deactivate_isolated_components!(data)
        num_iterations += 1
    end

    Memento.info(_LOGGER, "Network simplification reached in $(num_iterations) rounds.")

    # Returns whether or not the data has been modified.
    return revised
end


"""
Propagates inactive network node statuses to attached components (e.g., pipes) so that
system status values are consistent. Returns true if any component was modified.
"""
function propagate_topology_status!(data::Dict{String, <:Any})::Bool
    revised = false
    wm_data = get_wm_data(data)

    if _IM.ismultinetwork(wm_data)
        for wm_nw_data in values(wm_data["nw"])
            revised |= _propagate_topology_status!(wm_nw_data)
        end
    else
        revised = _propagate_topology_status!(wm_data)
    end

    return revised
end


""
function _propagate_topology_status!(data::Dict{String,<:Any})
    nodes = Dict{Int, Any}(node["index"] => node for (i, node) in data["node"])

    # Compute what active demands are incident to each node.
    incident_demand = node_comp_lookup(data["demand"], data["node"])
    incident_active_demand = Dict()

    for (i, demand_list) in incident_demand
        incident_active_demand[i] = [demand for demand in demand_list if demand["status"] != STATUS_INACTIVE]
    end

    # Compute what active reservoirs are incident to each node.
    incident_reservoir = node_comp_lookup(data["reservoir"], data["node"])
    incident_active_reservoir = Dict()

    for (i, reservoir_list) in incident_reservoir
        incident_active_reservoir[i] = [reservoir for reservoir in reservoir_list if reservoir["status"] != STATUS_INACTIVE]
    end

    # Compute what active tanks are incident to each node.
    incident_tank = node_comp_lookup(data["tank"], data["node"])
    incident_active_tank = Dict()

    for (i, tank_list) in incident_tank
        incident_active_tank[i] = [tank for tank in tank_list if tank["status"] != STATUS_INACTIVE]
    end

    # Compute what active pipes are incident to each node.
    incident_pipe = Dict{Int, Any}(node["index"] => [] for (i, node) in data["node"])

    for pipe in values(data["pipe"])
        push!(incident_pipe[pipe["node_fr"]], pipe)
        push!(incident_pipe[pipe["node_to"]], pipe)
    end

    # Compute what active design pipes are incident to each node.
    incident_des_pipe = Dict{Int, Any}(node["index"] => [] for (i, node) in data["node"])

    for des_pipe in values(data["des_pipe"])
        push!(incident_des_pipe[des_pipe["node_fr"]], des_pipe)
        push!(incident_des_pipe[des_pipe["node_to"]], des_pipe)
    end

    # Compute what active pumps are incident to each node.
    incident_pump = Dict{Int, Any}(node["index"] => [] for (i, node) in data["node"])

    for pump in values(data["pump"])
        push!(incident_pump[pump["node_fr"]], pump)
        push!(incident_pump[pump["node_to"]], pump)
    end

    # Compute what active expansion pumps are incident to each node.
    incident_ne_pump = Dict{Int, Any}(node["index"] => [] for (i, node) in data["node"])

    for ne_pump in values(data["ne_pump"])
        push!(incident_ne_pump[ne_pump["node_fr"]], ne_pump)
        push!(incident_ne_pump[ne_pump["node_to"]], ne_pump)
    end

    # Compute what active regulators are incident to each node.
    incident_regulator = Dict{Int, Any}(node["index"] => [] for (i, node) in data["node"])

    for regulator in values(data["regulator"])
        push!(incident_regulator[regulator["node_fr"]], regulator)
        push!(incident_regulator[regulator["node_to"]], regulator)
    end

    # Compute what active short pipes are incident to each node.
    incident_short_pipe = Dict{Int, Any}(node["index"] => [] for (i, node) in data["node"])

    for short_pipe in values(data["short_pipe"])
        push!(incident_short_pipe[short_pipe["node_fr"]], short_pipe)
        push!(incident_short_pipe[short_pipe["node_to"]], short_pipe)
    end

    # Compute what active network expansion short pipes are incident to each node.
    incident_ne_short_pipe = Dict{Int, Any}(node["index"] => [] for (i, node) in data["node"])

    for ne_short_pipe in values(data["ne_short_pipe"])
        push!(incident_ne_short_pipe[ne_short_pipe["node_fr"]], ne_short_pipe)
        push!(incident_ne_short_pipe[ne_short_pipe["node_to"]], ne_short_pipe)
    end

    # Compute what active valves are incident to each node.
    incident_valve = Dict{Int, Any}(node["index"] => [] for (i, node) in data["node"])

    for valve in values(data["valve"])
        push!(incident_valve[valve["node_fr"]], valve)
        push!(incident_valve[valve["node_to"]], valve)
    end

    revised = false

    for comp_type in _LINK_COMPONENTS
        for (i, comp) in data[comp_type]
            if comp["status"] != STATUS_INACTIVE
                node_fr = nodes[comp["node_fr"]]
                node_to = nodes[comp["node_to"]]

                if any(x["status"] == STATUS_INACTIVE for x in [node_fr, node_to])
                    message = "Deactivating $(replace(comp_type, "_" => " ")) $(i): "
                    message *= "($(comp["node_fr"]), $(comp["node_to"])) "
                    message *= "because of a connecting node's status."
                    Memento.info(_LOGGER, message)

                    comp["status"] = STATUS_INACTIVE
                    revised = true
                end
            end
        end
    end

    for (i, node) in nodes
        if node["status"] == STATUS_INACTIVE
            for demand in incident_active_demand[i]
                if demand["status"] != STATUS_INACTIVE
                    message = "Deactivating demand $(demand["index"]) due to inactive node $(i)."
                    Memento.info(_LOGGER, message)
                    demand["status"] = STATUS_INACTIVE
                    revised = true
                end
            end

            for reservoir in incident_active_reservoir[i]
                if reservoir["status"] != STATUS_INACTIVE
                    message = "Deactivating reservoir $(reservoir["index"]) due to inactive node $(i)."
                    Memento.info(_LOGGER, message)
                    reservoir["status"] = STATUS_INACTIVE
                    revised = true
                end
            end

            for tank in incident_active_tank[i]
                if tank["status"] != STATUS_INACTIVE
                    message = "Deactivating tank $(tank["index"]) due to inactive node $(i)."
                    Memento.info(_LOGGER, message)
                    tank["status"] = STATUS_INACTIVE
                    revised = true
                end
            end
        end
    end

    return revised
end


"Builds a lookup list of what components are connected to a given node."
function node_comp_lookup(comp_data::Dict{String,<:Any}, node_data::Dict{String,<:Any})
    node_comp = Dict{Int, Any}(node["index"] => [] for (i, node) in node_data)

    for comp in values(comp_data)
        push!(node_comp[comp["node"]], comp)
    end

    return node_comp
end
