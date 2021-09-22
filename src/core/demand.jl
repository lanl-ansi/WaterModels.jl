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


function make_demands_dispatchable!(data::Dict{String,<:Any})
    apply_wm!(_make_demands_dispatchable!, data)
end


function _make_demands_dispatchable!(data::Dict{String,<:Any})
    _make_demand_dispatchable!.(values(data["demand"]))
end


function _make_demand_dispatchable!(demand::Dict{String,<:Any})
    if demand["flow_min"] >= 0.0 && demand["flow_max"] > 0.0
        # If the demand is a sink for flow, modify the lower bound.
        demand["flow_min"] = 0.0
        demand["dispatchable"] = true
    elseif demand["flow_min"] < 0.0 && demand["flow_max"] <= 0.0
        # If the demand is a source of flow, modify the upper bound.
        demand["flow_max"] = 0.0
        demand["dispatchable"] = true
    end
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
    apply_wm!(_set_demand_warm_start!, data)
end


function _set_demand_warm_start!(data::Dict{String, <:Any})
    for demand in values(data["demand"])
        flow_mid = 0.5 * (demand["flow_min"] + demand["flow_max"])
        demand["q_demand_start"] = get(demand, "q", flow_mid)
    end
end
