# using LightGraphs
# using MetaGraphs
# using Metis

# using GraphPlot
# using Cairo, Compose

function build_graph_meta(data::Dict{String, <:Any})
    num_nodes = length(data["node"])
    graph = LightGraphs.SimpleGraph(num_nodes)
    Gr    = MetaGraph(graph, 100.0)
    node_names = sort(collect(keys(data["node"])))
    node_map = Dict{String,Int}(x => i for (i, x) in enumerate(node_names))

    weights = []
    for component_type in ["pipe", "pump", "regulator", "short_pipe", "valve"]
        for (a, comp) in data[component_type]
            i, j = node_map[string(comp["node_fr"])], node_map[string(comp["node_to"])]
            MetaGraphs.add_edge!(Gr, i, j)
            if component_type == "pipe"
                set_prop!(Gr, i, j, :weight, 1)
                append!(weights, 1)
            else
                set_prop!(Gr, i, j, :weight, 100)
                append!(weights, 100)
            end

        end
    end

    return Gr, node_map, weights
end

function build_graph(data::Dict{String, <:Any})
    num_nodes = length(data["node"])
    graph = LightGraphs.SimpleGraph(num_nodes)
    Gr    = Graph(graph)
    node_names = sort(collect(keys(data["node"])))
    node_map = Dict{String,Int}(x => i for (i, x) in enumerate(node_names))

    weights = []
    for component_type in ["pipe", "pump", "regulator", "short_pipe", "valve"]
        for (a, comp) in data[component_type]
            i, j = node_map[string(comp["node_fr"])], node_map[string(comp["node_to"])]
            LightGraphs.add_edge!(Gr, i, j)
            if component_type == "pipe"
                append!(weights, 1)
            else
                append!(weights, 100)
            end
        end
    end

    return Gr, node_map, weights
end

function create_node_partition(data::Dict{String, <:Any}, num_partitions)
    G, node_map, w = build_graph_meta(data)
    

end
