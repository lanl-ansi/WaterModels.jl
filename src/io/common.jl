"""
    parse_file(path)

Parses an [EPANET](https://www.epa.gov/water-research/epanet) (.inp) or
JavaScript Object Notation (JSON) file from the file path `path`, depending on
the file extension, and returns a WaterModels data structure (a dictionary of
data).
"""
function parse_file(path::String)
    if endswith(path, ".inp")
        network_data = WaterModels.parse_epanet(path)
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
    dict = JSON.parsefile(path)
    dict["per_unit"] = false
    return dict
end
