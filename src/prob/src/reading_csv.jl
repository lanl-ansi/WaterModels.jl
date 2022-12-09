max_out_junc_number = 1000
max_out_edge_number = 1000
max_pl_edge_number = max_out_edge_number*10
function parse_csv(junc_file::AbstractString,pipe_file::AbstractString, pl_file::AbstractString)
    open(junc_file, "r") do junc_io
        open(pipe_file, "r") do pipe_io
            open(pl_file, "r") do pl_io
                junc_data = parse_junc_csv(junc_io)
                pipe_data = parse_pipe_csv(pipe_io)
                (plant_data, load_data) = parse_plant_load_csv(pl_io)
                data = Dict()
                get!(data, "junctions", junc_data)
                get!(data, "pipes", pipe_data)
                get!(data, "plants", plant_data)
                get!(data, "loads", load_data)
                return data
            end
        end
    end
end

function parse_junc_csv(io::IO)::Array{Dict{String,Any},1}
    raw = readlines(io)

    data = []
    for line in raw[2:end]
        id,x_loc, y_loc = split(line, ",")
        push!(
            data,
            Dict(
                "id" => parse(Int64, id),
                "type" => 'o',
                "x_location" => x_loc,
                "y_location" => y_loc,
            ),
        )
        push!(
            data,
            Dict(
                "id" => parse(Int, id) + max_out_junc_number,
                "type" => 'r',
                "x_location" => x_loc,
                "y_location" => y_loc,
            ),
        )
    end
    return data
end

function parse_pipe_csv(edge_io::IO)::Array{Dict{String,Any},1}
    raw = readlines(edge_io)

    pipe_data = []
    for line in raw[2:end]
        id, fr_node, to_node, out_diameter, return_diameter, length = split(line, ",")
        push!(
            pipe_data,
            Dict(
                "id" => parse(Int, id),
                "type" => 'o',
                "fr_node" => parse(Int, fr_node),
                "to_node" => parse(Int, to_node),
                "diameter" => parse(Float64, out_diameter),
                "length" => parse(Float64, length),
            ),
        )
        push!(
            pipe_data,
            Dict(
                "id" => parse(Int, id) + max_out_edge_number,
                "type" => 'r',
                "fr_node" => parse(Int, fr_node) + max_out_junc_number,
                "to_node" => parse(Int, to_node) + max_out_junc_number,
                "diameter" => parse(Float64,return_diameter),
                "length" => parse(Float64, length),
            ),
        )
    end
    return pipe_data
end

function  parse_plant_load_csv(pl_io::IO)
    raw = readlines(pl_io)

    plant_data = [];
    load_data = [];
    for line in raw[2:end]
        pl_id, node_id, name, load = split(line, ",")
        id = parse(Int, pl_id) + max_pl_edge_number
        fr_node = parse(Int, node_id)
        to_node = parse(Int, node_id) + max_out_edge_number
        load = parse(Float64, load)

        if(name == "Plant")
            push!(
                plant_data,
                Dict(
                "id" => id,
                "fr_node" => fr_node,
                "to_node" => to_node,
                "load" => load,
                ),
            )
        else
            push!(
                load_data,
                Dict(
                "id" => id,
                "fr_node" => fr_node,
                "to_node" => to_node,
                "load" => load,
                ),
            )
        end
    end
    return plant_data, load_data
end
