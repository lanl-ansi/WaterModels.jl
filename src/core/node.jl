function _relax_nodes!(data::Dict{String,<:Any})
    if !_IM.ismultinetwork(data)
        ts_keys = keys(data["time_series"])

        if "node" in keys(data["time_series"]) && isa(data["time_series"]["node"], Dict)
            ts = data["time_series"]["node"]
            nodes = values(filter(x -> x.first in keys(ts), data["node"]))
            map(x -> x["h_min"] = minimum(ts[string(x["index"])]["head"]), nodes)
            map(x -> x["h_max"] = maximum(ts[string(x["index"])]["head"]), nodes)
        end
    end
end
