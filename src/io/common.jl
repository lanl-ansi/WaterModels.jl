"""
    parse_file(
        path::String;
        skip_correct::Bool=false,
        per_unit::Bool=true
    )

Parses an [EPANET](https://www.epa.gov/water-research/epanet) (.inp) or
JavaScript Object Notation (JSON) file from the file path `path`, depending
on the file extension, and returns a WaterModels data structure (i.e., a
dictionary of data). Here, `skip_correct` will skip data correction routines
(e.g., component status propagation) if set to `true`, and `per_unit` will
translate the data model to a per-unit measurement system if set to `true`.
"""
function parse_file(path::String; skip_correct::Bool = false, per_unit::Bool = true)
    if endswith(path, ".inp")
        network_data = WaterModels.parse_epanet(path)
    elseif endswith(path, ".json")
        network_data = WaterModels.parse_json(path)
        correct_enums!(network_data)
    else
        error("\"$(path)\" is not a valid file type.")
    end

    if !haskey(network_data, "per_unit")
        network_data["per_unit"] = false
    end

    if !skip_correct
        correct_network_data!(network_data; per_unit=per_unit)
    end

    return network_data
end


"""
    parse_json(path::String)

Parses a JavaScript Object Notation (JSON) file from the file path `path` and
returns a WaterModels data structure (i.e., a dictionary of data). Does not
perform data correction nor per-unit translations of the data model.
"""
function parse_json(path::String)
    return JSON.parsefile(path)
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


function correct_network_data!(data::Dict{String, <:Any}; per_unit::Bool = true)
    epanet_to_watermodels!(data; import_all = false)
    correct_pipes!(data)
    correct_des_pipes!(data)
    correct_pumps!(data)
    correct_regulators!(data)
    correct_short_pipes!(data)
    correct_valves!(data)
    correct_nodes!(data)

    if per_unit
        # Make data per-unit if necessary.
        make_per_unit!(data)
    end
end
