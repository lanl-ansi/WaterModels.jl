function _relax_demand!(demand::Dict{String,<:Any})
    if haskey(demand, "demand_min") && haskey(demand, "demand_max")
        demand["dispatchable"] = demand["demand_min"] != demand["demand_max"]
    else
        demand["dispatchable"] = false
    end
end


function _relax_demands!(data::Dict{String,<:Any})
    if _IM.ismultinetwork(data)
        for (n, nw) in data["nw"]
            map(x -> _relax_demand!(x), values(nw["demand"]))
        end
    else
        if haskey(data, "time_series") && haskey(data["time_series"], "demand")
            ts = data["time_series"]["demand"]
            dems = values(filter(x -> x.first in keys(ts), data["demand"]))
            map(x -> x["demand_min"] = minimum(ts[string(x["index"])]["flow_rate"]), dems)
            map(x -> x["demand_max"] = maximum(ts[string(x["index"])]["flow_rate"]), dems)
        end

        map(x -> _relax_demand!(x), values(data["demand"]))
    end
end
