# reading junction data
n_junctions = size(junctions_mat)[1];
junctions_dict = Dict();

for i in 1:n_junctions
    id = junctions_mat[i,1]
    junc_dict = Dict();
    get!(junc_dict, "id", junctions_mat[i,1])
    get!(junc_dict, "type", junctions_mat[i,2])
    get!(junctions_dict, id, junc_dict)
end

SetNO = Set(); SetNR = Set();
for (id,dict) in junctions_dict
    if(dict["type"] == 'o')
        union!(SetNO, dict["id"])
    elseif(dict["type"] == 'r')
        union!(SetNR, dict["id"])
    end
end

SetN = union(SetNO, SetNR);

#reading pipe data
n_pipes = size(pipes_mat)[1];
pipes_dict = Dict();

for i in 1:n_pipes
    id = pipes_mat[i,1]
    pipe_dict = Dict();
    get!(pipe_dict, "id", pipes_mat[i,1])
    get!(pipe_dict, "type", pipes_mat[i,2])
    get!(pipe_dict, "fr_junction", pipes_mat[i,3])
    get!(pipe_dict, "to_junction", pipes_mat[i,4])
    get!(pipe_dict, "diameter", pipes_mat[i,5])
    get!(pipe_dict, "length", pipes_mat[i,6])
    get!(pipe_dict, "friction_factor", pipes_mat[i,7])
    get!(pipe_dict, "thermal_loss_coefficient", pipes_mat[i,8])
    get!(pipes_dict, id, pipe_dict)
end

SetPO = Set(); SetPR = Set();
pipes_entering_junction = Dict();
pipes_leaving_junction = Dict();

for (id, dict) in junctions_dict
    get!(pipes_leaving_junction, id, Set())
    get!(pipes_entering_junction, id, Set())
end

for (id, dict) in pipes_dict
    if(haskey(pipes_entering_junction, dict["to_junction"]))
        pipes_entering = union!(pipes_entering_junction[dict["to_junction"]], dict["id"])
    else
        pipes_entering = dict["id"]
    end
    get!(pipes_entering_junction, dict["to_junction"], pipes_entering)
    if(haskey(pipes_leaving_junction, dict["fr_junction"]))
        pipes_leaving = union!(pipes_leaving_junction[dict["fr_junction"]], dict["id"])
    else
        pipes_leaving = dict["id"]
    end
    get!(pipes_leaving_junction, dict["fr_junction"], pipes_leaving)
    if(dict["type"] == 'o')
        union!(SetPO, dict["id"])
    elseif(dict["type"]=='r')
        union!(SetPR, dict["id"])
    end
end
SetP = union(SetPO, SetPR)

#reading plant data
n_plants = size(plants_mat)[1];
plants_dict = Dict();

for i in 1:n_plants
    id = plants_mat[i,1]
    plant_dict = Dict();
    get!(plant_dict, "id", Int(plants_mat[i,1]))
    get!(plant_dict, "fr_junction", Int(plants_mat[i,2]))
    get!(plant_dict, "to_junction", Int(plants_mat[i,3]))
    get!(plants_dict, id, plant_dict)
end

PlantSet = Set();
plants_entering_junction = Dict();
plants_leaving_junction = Dict();
for (id, dict) in junctions_dict
    get!(plants_leaving_junction, id, Set())
    get!(plants_entering_junction, id, Set())
end

for(id, dict) in plants_dict
    union!(PlantSet, dict["id"])
    if(haskey(plants_entering_junction, dict["to_junction"]))
        plants_entering = union!(plants_entering_junction[dict["to_junction"]], dict["id"])
    else
        plants_entering = dict["id"]
    end
    get!(plants_entering_junction, dict["to_junction"], dict["id"])

    if(haskey(plants_leaving_junction, dict["fr_junction"]))
        plants_leaving = union!(plants_leaving_junction[dict["fr_junction"]], dict["id"])
    else
        plants_leaving = dict["id"]
    end
    get!(plants_leaving_junction, dict["fr_junction"], dict["id"])
end

#reading load data
n_loads = size(loads_mat)[1];
loads_dict = Dict();

for i in 1:n_loads
    id = loads_mat[i,1]
    load_dict = Dict();
    get!(load_dict, "id", Int(loads_mat[i,1]))
    get!(load_dict , "fr_junction", Int(loads_mat[i, 2]))
    get!(load_dict, "to_junction", Int(loads_mat[i,3]))
    get!(load_dict, "Q", loads_mat[i,4])
    get!(loads_dict, id, load_dict)
end

LoadSet = Set();
loads_entering_junction = Dict();
loads_leaving_junction = Dict();
for (id, dict) in junctions_dict
    get!(loads_leaving_junction, id, Set())
    get!(loads_entering_junction, id, Set())
end

for (id, dict) in loads_dict
    union!(LoadSet, dict["id"])
    if(haskey(loads_entering_junction, dict["to_junction"]))
        loads_entering = union!(loads_entering_junction[dict["to_junction"]], dict["id"])
    else
        loads_entering = dict["id"]
    end
    get!(loads_entering_junction, dict["to_junction"], loads_entering)
    if(haskey(loads_leaving_junction, dict["to_junction"]))
        loads_leaving = union!(loads_leaving_junction[dict["fr_junction"]], dict["id"])
    else
        loads_leaving = dict["id"]
    end
    get!(loads_leaving_junction, dict["fr_junction"], loads_leaving)
end

#Combined Edges
n_edges = n_pipes+n_plants+n_loads
SetE = union(SetPO, SetPR, LoadSet, PlantSet)
edges_dict = merge(pipes_dict, plants_dict, loads_dict)


edges_entering_junction = Dict();
edges_leaving_junction = Dict();
# for (id, dict) in junctions_dict
#     get!(edges_leaving_junction, id, Set())
#     get!(edges_entering_junction, id, Set())
# end

for (id, dict) in junctions_dict
    pipes_entering = Set(); plants_entering = Set(); loads_entering = Set();
    if(haskey(pipes_entering_junction, id))
        pipes_entering = pipes_entering_junction[id]
    end
    if(haskey(plants_entering_junction, id))
        plants_entering = plants_entering_junction[id]
    end
    if(haskey(loads_entering_junction, id))
        loads_entering = loads_entering_junction[id]
    end
    edges_entering = union(pipes_entering, plants_entering, loads_entering)
    get!(edges_entering_junction, id, edges_entering)

    pipes_leaving = Set(); plants_leaving = Set(); loads_leaving = Set();
    pipes_leaving = Set(); plants_entering = Set(); loads_entering = Set();
    if(haskey(pipes_leaving_junction, id))
        pipes_leaving = pipes_leaving_junction[id]
    end
    if(haskey(plants_leaving_junction, id))
        plants_leaving = plants_leaving_junction[id]
    end
    if(haskey(loads_leaving_junction, id))
        loads_leaving = loads_leaving_junction[id]
    end
    edges_leaving = union(pipes_leaving, plants_leaving, loads_leaving)
    get!(edges_leaving_junction, id, edges_leaving)
end

################################################################################
E_in = edges_entering_junction
E_out = edges_leaving_junction

P_in = pipes_entering_junction
P_out = pipes_leaving_junction

PL_in = plants_entering_junction
PL_out = plants_leaving_junction

L_in = loads_entering_junction
L_out = loads_leaving_junction


# N_to = zeros(Int,n_edges)
# N_fr = zeros(Int,n_edges)
N_to = Dict()
N_fr = Dict()
for e in SetE
    get!(N_to,e,edges_dict[e]["to_junction"])
    get!(N_fr,e,edges_dict[e]["fr_junction"])
    # N_to[e] = edges_dict[e]["to_junction"]
    # N_fr[e] = edges_dict[e]["fr_junction"]
end

plant_froms = Set()
for (id,dict) in plants_dict
    push!(plant_froms,dict["fr_junction"])
end
