function _relax_reservoirs!(data::Dict{String,<:Any})
    if _IM.ismultinetwork(data)
        reservoirs = vcat([vcat(values(nw["reservoir"])...) for (n, nw) in data["nw"]]...)
        map(x -> x["dispatchable"] = true, reservoirs)
    else
        reservoirs = values(data["reservoir"])
        map(x -> x["dispatchable"] = true, reservoirs)
    end
end
