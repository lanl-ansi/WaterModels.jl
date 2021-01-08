"""
    parse_file(path)

Parses an [EPANET](https://www.epa.gov/water-research/epanet) (.inp) or JavaScript Object
Notation (JSON) file from the file path `path`, depending on the file extension, and returns
a WaterModels data structure (a dictionary of data).
"""
function parse_file(path::String; skip_correct::Bool = false)
    if endswith(path, ".inp")
        network_data = WaterModels.parse_epanet(path)
        epanet_to_watermodels!(network_data; import_all = false)
    else
        error("\"$(path)\" is not a valid file type.")
    end

    if !skip_correct
        correct_network_data!(network_data)
    end

    return network_data
end


"""
    parse_json(path)

Parses a JavaScript Object Notation (JSON) file from the file path `path` and returns a
WaterModels data structure (a dictionary of data).
"""
function parse_json(path::String)
    dict = JSON.parsefile(path)

    if !haskey(dict, "per_unit")
        dict["per_unit"] = false
    end

    return dict
end


"""
    _read_file_as_string(path)
"""
function _read_file_as_string(file_path::String)
    if isfile(file_path)
        return read(open(file_path), String)
    else
        error_message = "File \"$(file_path)\" does not exist."
        Memento.error(_LOGGER, error_message)
    end
end


function correct_network_data!(data::Dict{String, <:Any})
    # Correct edge-type component data.
    correct_pipes!(data)
    correct_des_pipes!(data)
    correct_pumps!(data)
    correct_regulators!(data)
    correct_short_pipes!(data)
    correct_valves!(data)

    # Correct nodal component data.
    correct_nodes!(data)
end
