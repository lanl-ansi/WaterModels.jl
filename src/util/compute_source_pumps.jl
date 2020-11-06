"""
This function takes in the watermodels data dictionary and returns a dictionary of all directional links in
the water network that result from components on the network.

"""
function find_directional_components(data::Dict{String,<:Any})

    directional_components = Dict()
    component_id  = 0 

    for component_type in ["pipe", "pump", "regulator", "short_pipe", "valve"]
        if component_type == "regulator"
            for (a, comp) in data[component_type]
                component_id += 1
                i, j = comp["node_fr"], comp["node_to"]
                directional_components[component_id] = [i,j]
            end
        else
            for (a, comp) in data[component_type]
                if comp["flow_direction"] == POSITIVE
                    component_id += 1
                    #i, j = node_map[string(comp["node_fr"])], node_map[string(comp["node_to"])]
                    i, j = comp["node_fr"], comp["node_to"]
                    directional_components[component_id] = [i,j]
                elseif  comp["flow_direction"] == NEGATIVE
                    i, j = comp["node_fr"], comp["node_to"]
                    directional_components[component_id] = [j,i]

                
                end

            end
        end
    end
    return directional_components
end 

"""
This function takes in the watermodels data dictionary and returns a dictionary of all directional links associated 
with pumps in the water network.

"""
function find_directional_pumps(data::Dict{String,<:Any})
    pumps_dict = Dict()
    #pump_id = 0
  
    for (a, comp) in data["pump"]
        if comp["flow_direction"] == POSITIVE
            #pump_id += 1
            #i, j = node_map[string(comp["node_fr"])], node_map[string(comp["node_to"])]
            i, j = comp["node_fr"], comp["node_to"]
            pumps_dict[a] = [i,j]
        elseif  comp["flow_direction"] == NEGATIVE
            i, j = comp["node_fr"], comp["node_to"]
            pumps_dict[a] = [j,i]
        end
    end
    return pumps_dict
end 

"""
This function checks if the edge from current node to next node if deleted disconnects the graph. 
If the graph is disconnected it then checks if next_node is in the same component of the graph as sink.
If not the edge fails the test. Otherwise the edge passes.
"""
function component_check(sink::Int,current_node::Int,next_node::Int,graph::LightGraphs.SimpleGraphs.AbstractSimpleGraph{Int64})
    local current_node_comp
    local sink_node_comp
    test_graph = copy(graph)

    LightGraphs.rem_edge!(test_graph, current_node,next_node)

    cc_graph = LightGraphs.connected_components(graph)
    cc_test_graph = LightGraphs.connected_components(test_graph)

    if length(cc_graph) < length(cc_test_graph) 
        for comp in cc_test_graph
            if next_node in comp
                current_node_comp = comp
                break
            end  
        end
        for comp in cc_test_graph
            if sink in comp
                sink_node_comp = comp
                break
            end  
        end

        if sink_node_comp  == current_node_comp
            return "passed"
        else
            return "failed"
        end
    end

    return "passed"
end
"""
This function finds all paths from the start node to the end node in a graph.
"""
function find_paths(start::Int, _end::Int, graph::LightGraphs.SimpleGraphs.AbstractSimpleGraph{Int64})

    #intialize data structures
    paths_found = 0
    paths_dict = Dict()
    path = [] 

    # dictionary for tracking remaining possible movements
    poss = Dict()

    for i in range(1,length=length(graph.fadjlist))
        poss[i] = []
        for j in graph.fadjlist[i]
            push!(poss[i],j) 
        end
    end

    # start algorithm
    push!(path,start);

    while isempty(path) == false

        tail = last(path)

        if isempty(poss[tail]) == false

            next = pop!(poss[tail]) # go to next and remove it as a movement option from tail
       
            if next == _end
                push!(path,next)
                paths_found += 1
                paths_dict[paths_found] = copy(path)

                # set up for next attempt
                pop!(path)
                filter!(x ->(x!=_end),poss[last(path)]) 
          
            elseif next in path
                filter!(x ->(x!=next),poss[last(path)])

            elseif component_check(_end,tail,next,graph) == "failed"
                filter!(x ->(x!=next),poss[last(path)])

            else
                push!(path,next) # next vertex onto path
                filter!(x -> (x!=tail),poss[next])

            end
 
        else
            v = pop!(path)
            if isempty(path) == true
                break
            end

            filter!(x ->(x!=v),poss[last(path)])
            poss[v] = copy(graph.fadjlist[v])
   
        end

    end

    return paths_dict
end
"""
This function checks to see if a path in the graph has a pump on it. 
"""
function check_path_for_pump(path::Array{Any,1},pumps_dict::Dict{Any,Any},node_map::Dict{String,Int64})
    status = false
    for k in range(1, stop = length(path) - 1 )
        i,j = path[k],path[k+1]
        for (pump_id,pump) in pumps_dict

            pump_start,pump_end = node_map[string(pump[1])],node_map[string(pump[2])]

            if [pump_start,pump_end] == [i,j]
                status = true
                break
            end
        end
        if status == true
            break
        end
    
    end 

    return status
end

"""
This function checks to see if a path in the graph has a directional component on it making it so water
could not flow from the start of the path to end of the path along the path. 
"""
function check_path_for_directional_components(path::Array{Any,1},directional_components_dict::Dict{Any,Any},node_map::Dict{String,Int64})
    status = false
    for k in range(1, stop = length(path) - 1)
        i,j = path[k],path[k+1]
        for (component_id,component) in directional_components_dict

            component_start,component_end = node_map[string(component[1])],node_map[string(component[2])]

            if [component_start,component_end] == [j,i]
                status = true
                break
            end
        end
        if status == true
            break
        end
    
    end 

    return status
end

# function check_elevation_difference(path,network,node_map)
#     status = "downhill"
#     local start_path 
#     local end_path 
#     for (key,val) in node_map

#         if val == last(path)
#             end_path = key
#         end

#         if val == first(path)
#             start_path = key
#         end
#     end
    
#     start_elevation = network["node"][start_path]["elevation"]
#     end_elevation = network["node"][end_path]["elevation"]

#     if end_elevation - start_elevation > 0 
#         status = "uphill"
#     end
    
#     return status,end_elevation - start_elevation
# end
"""
This function checks to see if a path from a reservoir to the start of a pump in the graph has a directional component or pumps on it making it so water
flowing from the start of the path to end of the path along the path could either not happen or would be water that has already passed through a pump. 
"""
function check_path(path::Array{Any,1},pumps_dict::Dict{Any,Any},directional_components_dict::Dict{Any,Any},node_map::Dict{String,Int64})

    
    pump_status = check_path_for_pump(path,pumps_dict,node_map)
    if pump_status == true
        return "has pump"
    end

    dc_status = check_path_for_directional_components(path,directional_components_dict,node_map)
    if dc_status == true
        return "not a path"
    end

    # elevation_status, diff = check_elevation_difference(path,network,node_map)
    # if elevation_status == "uphill"
    #     return "uphill path", diff
    # end

    return "valid path"

end

"""
This function checks to see if a pump when on only pulls water from reservoirs. 
"""
function check_if_source_pump(pump::Array{Int64,1},node_map::Dict{String,Int64},network::Dict{String,
    <:Any},pumps_dict::Dict{Any,Any},directional_components_dict::Dict{Any,Any},graph::LightGraphs.SimpleGraphs.AbstractSimpleGraph{Int64})
    pump_from_node = node_map[string(first(pump))]
    for (res_id,res) in network["reservoir"]
        res_node = node_map[string(res["node"])]
        paths_dict = find_paths(res_node, pump_from_node, graph)
        possible_paths = collect(keys(paths_dict))
        for (path_id,path) in paths_dict
            status = check_path(path,pumps_dict,directional_components_dict,node_map)
            if status == "valid path"
                return true 
            end
        end
    end
    return false
end
"""
This function goes through all the pumps on the water network and determines which ones are source pumps and which ones are not
and returns a list of source pumps.
"""
function find_source_pumps(network::Dict{String,<:Any})

    graph, node_map = create_graph(network)
    pumps_dict = find_directional_pumps(network);
    directional_components_dict = find_directional_components(network) 
    
    source_pump_ids = []

    for (pump_id,pump) in pumps_dict
        status = check_if_source_pump(pump,node_map,network,pumps_dict,directional_components_dict,graph)
        if status == true
            push!(source_pump_ids,pump_id)
        end
    end
    return source_pump_ids
end   



