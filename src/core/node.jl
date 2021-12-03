function relax_nodes!(data::Dict{String,Any})
    apply_wm!(_relax_nodes!, data; apply_to_subnetworks = true)
end


function _relax_nodes!(data::Dict{String,<:Any})
    if !_IM.ismultinetwork(data)
        if haskey(data, "time_series") && haskey(data["time_series"], "node")
            ts = data["time_series"]["node"]
            nodes = values(filter(x -> x.first in keys(ts), data["node"]))
            map(x -> x["head_min"] = minimum(ts[string(x["index"])]["head_min"]), nodes)
            map(x -> x["head_max"] = maximum(ts[string(x["index"])]["head_max"]), nodes)
        end
    end
end


function _set_node_bounds_from_time_series!(node::Dict{String,<:Any}, time_series::Dict{String,<:Any})
    # Get the index of the node.
    node_index = string(node["index"])

    if haskey(time_series, "node") && haskey(time_series["node"], node_index)
        # Get the time series data corresponding to the node.
        node_time_series = time_series["node"][node_index]
        
        if haskey(node_time_series, "head_min")
            # Set minimum head to the minimum across all time.
            node["head_min"] = minimum(node_time_series["head_min"])
        end

        if haskey(node_time_series, "head_max")
            # Set the maximum head to the maximum across all time.
            node["head_max"] = maximum(node_time_series["head_max"])
        end

        if haskey(node_time_series, "head_nominal")
            # Set the nominal head to the mean across all time.
            node["head_nominal"] = Statistics.mean(node_time_series["head_nominal"])
        end

        # Ensure node values are bounded as expected.
        @assert node["head_min"] <= node["head_nominal"]
        @assert node["head_nominal"] <= node["head_max"]
    end
end


function correct_nodes!(data::Dict{String, <:Any})
    head_offset = _calc_head_offset(data)
    apply_wm!(x -> _correct_nodes!(x, head_offset), data; apply_to_subnetworks = true)
end


function _correct_nodes!(data::Dict{String, <:Any}, head_offset::Float64)
    # Compute a global estimate for maximum head.
    head_max = _calc_head_max(data)    

    for (idx, node) in data["node"]
        _correct_status!(node)
        demands = filter(x -> x.second["node"] == node["index"], data["demand"])
        reservoirs = filter(x -> x.second["node"] == node["index"], data["reservoir"])
        tanks = filter(x -> x.second["node"] == node["index"], data["tank"])
        _correct_node_head_bounds!(node, demands, reservoirs, tanks, head_offset, head_max)
    end
end


function _correct_node_head_bounds!(
    node::Dict{String, <:Any}, demands::Dict{String, <:Any},
    reservoirs::Dict{String, <:Any}, tanks::Dict{String, <:Any},
    head_offset::Float64, head_max::Float64)
    # Compute minimum and maximum head bounds for the node.
    node["head_min"] = _calc_node_head_min(node, demands, reservoirs, tanks, head_offset)
    node["head_max"] = _calc_node_head_max(node, demands, reservoirs, tanks, head_max)
    @assert node["head_min"] <= node["head_max"]
end


function _calc_node_head_min(
    node::Dict{String, <:Any}, demands::Dict{String, <:Any},
    reservoirs::Dict{String, <:Any}, tanks::Dict{String, <:Any},
    head_offset::Float64)
    # Get possible stored node head bound data.
    head_nominal = get(node, "head_nominal", -Inf)
    head_max = get(node, "head_max", -Inf)
    head_min_base = min(head_nominal, head_max, get(node, "head_min", -Inf))

    if length(tanks) > 0
        # Compute properties of the minimum tank level across connected tanks.
        min_level_min = minimum(x["min_level"] for (i, x) in tanks)
        min_level_max = maximum(x["min_level"] for (i, x) in tanks)

        # All minimum levels for tanks at the node should be equal.
        @assert min_level_min == min_level_max

        # Return the head associated with the minimum level.
        return max(node["elevation"] + min_level_min, head_min_base)
    elseif length(tanks) + length(demands) + length(reservoirs) == 0
        return max(node["elevation"] + head_offset, head_min_base)
    else
        return max(node["elevation"], head_min_base)
    end
end


function _calc_node_head_max(
    node::Dict{String, <:Any}, demands::Dict{String, <:Any},
    reservoirs::Dict{String, <:Any}, tanks::Dict{String, <:Any}, head_max::Float64)
    # Get possible stored node head bound data.
    head_nominal = get(node, "head_nominal", Inf)
    head_min = get(node, "head_min", Inf)
    head_max_base = max(head_nominal, head_min, get(node, "head_max", Inf))

    if length(tanks) > 0
        # Compute properties of the maximum tank level across connected tanks.
        max_level_min = minimum(x["max_level"] for (i, x) in tanks)
        max_level_max = maximum(x["max_level"] for (i, x) in tanks)

        # All maximum levels for tanks at the node should be equal.
        @assert max_level_min == max_level_max

        # Return the head associated with the minimum level.
        return min(node["elevation"] + max_level_min, head_max, head_max_base)
    elseif length(reservoirs) > 0
        if haskey(node, "head_max")
            return node["head_max"]
        else
            @assert haskey(node, "head_nominal")
            return get(node, "head_nominal", Inf)
        end
    else
        return min(head_max, head_max_base)
    end
end


function _calc_node_head_midpoint(node::Dict{String,Any})
    return 0.5 * (node["head_max"] + node["head_min"])
end


function calc_node_head_median_midpoint(data::Dict{String,Any})
    apply_wm!(_calc_node_head_median_midpoint, data; apply_to_subnetworks = true)
end


function _calc_node_head_median_midpoint(data::Dict{String,Any})
    if _IM.ismultinetwork(data)
        head_medians = Vector{Float64}([])

        for (n, nw) in data["nw"]
            head_midpoints = [_calc_node_head_midpoint(x) for (i, x) in nw["node"]]
            push!(head_medians, Statistics.median(head_midpoints)) # Add midpoint.
        end

        return Statistics.median(head_medians)
    else
        head_midpoints = [_calc_node_head_midpoint(x) for (i, x) in data["node"]]
        return Statistics.median(head_midpoints) # Calculate median midpoint.
    end
end


function set_node_warm_start!(data::Dict{String, <:Any})
    apply_wm!(_set_node_warm_start!, data)
end


function _set_node_warm_start!(data::Dict{String, <:Any})
    for node in values(data["node"])
        head_mid = 0.5 * (node["head_min"] + node["head_max"])
        node["h_start"] = get(node, "h", head_mid)
        node["p_start"] = get(node, "p", head_mid - node["elevation"])
    end
end
