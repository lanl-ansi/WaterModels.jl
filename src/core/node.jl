function _relax_nodes!(data::Dict{String,<:Any})
    if !_IM.ismultinetwork(data)
        if haskey(data, "time_series") && haskey(data["time_series"], "node")
            ts = data["time_series"]["node"]
            nodes = values(filter(x -> x.first in keys(ts), data["node"]))
            map(x -> x["h_min"] = minimum(ts[string(x["index"])]["head"]), nodes)
            map(x -> x["h_max"] = maximum(ts[string(x["index"])]["head"]), nodes)
        end
    end
end
