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


function correct_nodes!(data::Dict{String, <:Any})
    # Compute a global estimate for maximum head.
    head_max = _calc_head_max(data)

    for (idx, node) in data["node"]
        demands = filter(x -> x.second["node"] == node["index"], data["demand"])
        reservoirs = filter(x -> x.second["node"] == node["index"], data["reservoir"])
        tanks = filter(x -> x.second["node"] == node["index"], data["tank"])
        _correct_node_head_bounds!(node, demands, reservoirs, tanks, head_max)
    end
end


function _correct_node_head_bounds!(
    node::Dict{String, <:Any}, demands::Dict{String, <:Any},
    reservoirs::Dict{String, <:Any}, tanks::Dict{String, <:Any}, head_max::Float64)
    # Compute minimum and maximum head bounds for the node.
    node["head_min"] = _calc_node_head_min(node, demands, reservoirs, tanks)
    node["head_max"] = _calc_node_head_max(node, demands, reservoirs, tanks, head_max)
end


function _calc_node_head_min(
    node::Dict{String, <:Any}, demands::Dict{String, <:Any},
    reservoirs::Dict{String, <:Any}, tanks::Dict{String, <:Any})
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
        @assert haskey(node, "head_nominal")
        return get(node, "head_nominal", Inf)
    else
        return min(head_max, head_max_base)
    end
end
