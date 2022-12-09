# data_org_functions

EdgeInfo = []

for e in SetE
    e_dict = Dict(
    "edge_id" => e,
    "fr_junction" => fr[e],
    "to_junction" => to[e]
    )
    push!(EdgeInfo, e_dict)
end

function E_in(e)
    return EdgeInfo[e]
end
