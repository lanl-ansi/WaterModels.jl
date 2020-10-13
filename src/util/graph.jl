function create_graph(data::Dict{String,<:Any})
    num_nodes = length(data["node"])
    graph = LightGraphs.SimpleGraph(num_nodes)
    node_names = sort(collect(keys(data["node"])))
    node_map = Dict{String,Int}(x => i for (i, x) in enumerate(node_names))

    for component_type in ["pipe", "pump", "regulator", "short_pipe", "valve"]
        for (a, comp) in data[component_type]
            i, j = node_map[string(comp["node_fr"])], node_map[string(comp["node_to"])]
            LightGraphs.add_edge!(graph, i, j)
        end
    end

    return graph, node_map
end
