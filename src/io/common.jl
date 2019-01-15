import JSON

"""
    parse_file(path)

Parses an EPANET (.inp) or JavaScript Object Notation (JSON) file from the file
path `path` and returns a WaterModels data structure (a dictionary of data).
"""
function parse_file(path::String)
    if endswith(path, ".inp")
        network_data = WaterModels.parse_epanet_file(path)
    elseif endswith(path, ".json")
        network_data = parse_json(path)
    else
        error("\"$(path)\" is not a valid file type.")
    end
    
    return network_data
end

"""
    parse_json(path)

Parses a JavaScript Object Notation (JSON) file from the file path `path` and
returns a WaterModels data structure (a dictionary of data).
"""
function parse_json(path::String)
    return JSON.parsefile(path)
end

"""
    print_solution(solution)

Pretty-prints network optimization problem `solution` data.
"""
function print_solution(solution::Dict{String, Any})
    InfrastructureModels.print_summary(solution)
end
