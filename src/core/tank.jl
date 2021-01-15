function aggregate_tanks(subnetworks::Array{Dict{String, Any}, 1})
    tanks = deepcopy(subnetworks[1]["tank"])

    for (i, x) in tanks
        x["min_level"] = max_subnetwork_values(subnetworks, "tank", i, "min_level")
        x["max_level"] = min_subnetwork_values(subnetworks, "tank", i, "max_level")
        x["min_vol"] = max_subnetwork_values(subnetworks, "tank", i, "min_vol")
        x["init_level"] = min_subnetwork_values(subnetworks, "tank", i, "init_level")
        x["dispatchable"] = all_subnetwork_values(subnetworks, "tank", i, "dispatchable")
        x["status"] = any_subnetwork_values(subnetworks, "tank", i, "status")
    end

    return tanks
end


function _relax_tanks!(data::Dict{String,<:Any})
    if _IM.ismultinetwork(data)
        tanks = vcat([vcat(values(nw["tank"])...) for (n, nw) in data["nw"]]...)
        map(x -> x["dispatchable"] = true, tanks)
    else
        tanks = values(data["tank"])
        map(x -> x["dispatchable"] = true, tanks)
    end
end


function make_tank_start_dispatchable!(data::Dict{String,<:Any})
    if _IM.ismultinetwork(data)
        nw_ids = sort(collect(keys(data["nw"])))
        start_nw = string(sort([parse(Int, i) for i in nw_ids])[1])

        for (i, tank) in data["nw"][start_nw]["tank"]
            tank["dispatchable"] = true
        end
    else
        for (i, tank) in data["tank"]
            tank["dispatchable"] = true
        end
    end
end


function _fix_tank!(tank::Dict{String,<:Any})
    tank["dispatchable"] = false
end


function _fix_tanks!(data::Dict{String, <:Any})
    if _IM.ismultinetwork(data)
        for (n, nw) in data["nw"]
            map(x -> _fix_tank!(x), values(nw["tank"]))
        end
    else
        map(x -> _fix_tank!(x), values(data["tank"]))
    end
end