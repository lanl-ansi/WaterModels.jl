function aggregate_reservoirs(subnetworks::Array{Dict{String, Any}, 1})
    reservoirs = deepcopy(subnetworks[1]["reservoir"])

    for (i, x) in reservoirs
        x["head_nominal"] = min_subnetwork_values(subnetworks, "reservoir", i, "head_nominal")
        x["dispatchable"] = all_subnetwork_values(subnetworks, "reservoir", i, "dispatchable")
        x["status"] = any_subnetwork_values(subnetworks, "reservoir", i, "status")
    end

    return reservoirs
end


function _relax_reservoirs!(data::Dict{String,<:Any})
    if _IM.ismultinetwork(data)
        reservoirs = vcat([vcat(values(nw["reservoir"])...) for (n, nw) in data["nw"]]...)
        map(x -> x["dispatchable"] = true, reservoirs)
    else
        reservoirs = values(data["reservoir"])
        map(x -> x["dispatchable"] = true, reservoirs)
    end
end


function _fix_reservoir!(reservoir::Dict{String,<:Any})
    reservoir["dispatchable"] = false
end


function _fix_reservoirs!(data::Dict{String, <:Any})
    if _IM.ismultinetwork(data)
        for (n, nw) in data["nw"]
            map(x -> _fix_reservoir!(x), values(nw["reservoir"]))
        end
    else
        map(x -> _fix_reservoir!(x), values(data["reservoir"]))
    end
end


function set_reservoir_warm_start!(data::Dict{String, <:Any})
    InfrastructureModels.apply!(_set_reservoir_warm_start!, data, wm_it_name)
end


function _set_reservoir_warm_start!(data::Dict{String, <:Any})
    for reservoir in values(data["reservoir"])
        reservoir["q_reservoir_start"] = get(reservoir, "q", 0.0)
    end
end
