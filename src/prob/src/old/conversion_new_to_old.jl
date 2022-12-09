E_in = edges_entering_junction
E_out = edges_leaving_junction

P_in = pipes_entering_junction
P_out = pipes_leaving_junction

PL_in = plants_entering_junction
PL_out = plants_leaving_junction

L_in = loads_entering_junction
L_out = loads_leaving_junction


N_to = zeros(Int,n_edges)
N_fr = zeros(Int,n_edges)
for e in SetE
    N_to[e] = edges_dict[e]["to_junction"]
    N_fr[e] = edges_dict[e]["fr_junction"]
end

plant_froms = Set()
for (id,dict) in plants_dict
    push!(plant_froms,dict["fr_junction"])
end
