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
    end
end
