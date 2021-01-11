# Functions for working with WaterModels data elements.


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
    network_ids = sort([parse(Int, x) for x in keys(data["nw"])]; rev = true)
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
        head_gain = calc_pump_head_gain_max(pump, node_fr, node_to)
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
    return InfrastructureModels.make_multinetwork(data, union(global_keys, _wm_global_keys))
end


function split_multinetwork(data::Dict{String, <:Any}, nw_ids::Array{Array{String, 1}, 1})
    # Ensure the data has the multinetwork attribute.
    @assert _IM.ismultinetwork(data) == true

    # Get sub-multinetwork datasets indexed by the "nw" key.
    sub_mn = [Dict{String, Any}(i => data["nw"][i] for i in ids) for ids in nw_ids]

    # Get all data not associated with multinetwork components
    g_data = filter(x -> x.first != "nw", data)

    # Return the new, split multinetwork data dictionaries.
    return [merge(g_data, Dict{String, Any}("nw" => sub_mn[i])) for i in 1:length(nw_ids)]
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


function fix_all_indicators!(data::Dict{String,<:Any})
    set_start!(data, "pump", "status", "z_min")
    set_start!(data, "pump", "status", "z_max")
    set_start!(data, "regulator", "status", "z_min")
    set_start!(data, "regulator", "status", "z_max")
    set_start!(data, "valve", "status", "z_min")
    set_start!(data, "valve", "status", "z_max")
end


"Translate a multinetwork dataset to a snapshot dataset with dispatchable components."
function _relax_network!(data::Dict{String,<:Any})
    _relax_nodes!(data)
    _relax_tanks!(data)
    _relax_reservoirs!(data)
    _relax_demands!(data)
end
