function aggregate_demands(subnetworks::Array{Dict{String, Any}, 1})
    demands = deepcopy(subnetworks[1]["demand"])

    for (i, x) in demands
        x["flow_min"] = sum_subnetwork_values(subnetworks, "demand", i, "flow_min")
        x["flow_max"] = sum_subnetwork_values(subnetworks, "demand", i, "flow_max")
        x["flow_nominal"] = sum_subnetwork_values(subnetworks, "demand", i, "flow_nominal")
        x["dispatchable"] = all_subnetwork_values(subnetworks, "demand", i, "dispatchable")
        x["status"] = any_subnetwork_values(subnetworks, "demand", i, "status")
    end

    return demands
end


function _relax_demand!(demand::Dict{String,<:Any})
    if haskey(demand, "flow_min") && haskey(demand, "flow_max")
        demand["dispatchable"] = demand["flow_min"] != demand["flow_max"]
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
            map(x -> x["flow_min"] = minimum(ts[string(x["index"])]["flow_min"]), dems)
            map(x -> x["flow_max"] = maximum(ts[string(x["index"])]["flow_max"]), dems)
        end

        map(x -> _relax_demand!(x), values(data["demand"]))
    end
end

function _fix_demand!(demand::Dict{String,<:Any})
    demand["dispatchable"] = false
end

function _fix_demands!(data::Dict{String, <:Any})
    if _IM.ismultinetwork(data)
        for (n, nw) in data["nw"]
            map(x -> _fix_demand!(x), values(nw["demand"]))
        end
    else
        map(x -> _fix_demand!(x), values(data["demand"]))
    end
end


function set_demand_warm_start!(data::Dict{String, <:Any})
    InfrastructureModels.apply!(_set_demand_warm_start!, data, wm_it_name)
end


function _set_demand_warm_start!(data::Dict{String, <:Any})
    for demand in values(data["demand"])
        flow_mid = 0.5 * (demand["flow_min"] + demand["flow_max"])
        demand["q_demand_start"] = get(demand, "q", flow_mid)
    end
end
