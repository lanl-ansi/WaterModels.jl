function export_geojson(data::Dict{String, <:Any}, path::String)
    dict = export_geojson(data)

    open(path, "w") do f
        JSON.print(f, dict, 4)
    end
end

function export_geojson(data::Dict{String, <:Any})
    nodes = merge(data["junctions"], data["reservoirs"])
    node_features = build_nodes_geojson(nodes)
    links = merge(data["pipes"], data["valves"])
    link_features = build_links_geojson(links, nodes)
    features = vcat(node_features, link_features)
    return Dict("type" => "FeatureCollection", "features" => features)
end

function build_links_geojson(data::Dict{String, <:Any}, nodes::Dict{String, <:Any})
    link_array = []

    for (key, item) in data
        linestring = build_linestring_geojson(item, nodes)
        link_array = vcat(link_array, linestring)
    end

    return link_array
end

function build_nodes_geojson(data::Dict{String, <:Any})
    node_array = []

    for (key, item) in data
        point = build_point_geojson(item)
        node_array = vcat(node_array, point)
    end

    return node_array
end

function build_linestring_geojson(data::Dict{String, <:Any}, nodes::Dict{String, <:Any})
    geometry = Dict{String, Any}("type" => "LineString")

    if data["q"] >= 0.0
        from_x = nodes[data["node1"]]["x"]
        from_y = nodes[data["node1"]]["y"]
        to_x = nodes[data["node2"]]["x"]
        to_y = nodes[data["node2"]]["y"]
        geometry["coordinates"] = [[from_x, from_y], [to_x, to_y]]
    else
        from_x = nodes[data["node2"]]["x"]
        from_y = nodes[data["node2"]]["y"]
        to_x = nodes[data["node1"]]["x"]
        to_y = nodes[data["node1"]]["y"]
        geometry["coordinates"] = [[from_x, from_y], [to_x, to_y]]
    end

    properties = deepcopy(data)
    return Dict{String, Any}("type" => "Feature", "geometry" => geometry, "properties" => properties)
end

function build_point_geojson(data::Dict{String, <:Any})
    geometry = Dict{String, Any}("type" => "Point", "coordinates" => [data["x"], data["y"]])
    properties = deepcopy(data)
    pop!(properties, "x")
    pop!(properties, "y")
    return Dict{String, Any}("type" => "Feature", "geometry" => geometry, "properties" => properties)
end
