############################################################################
#                                                                          #
# This file provides functionality for interfacing with EPANET data files. #
# See https://github.com/OpenWaterAnalytics/EPANET/wiki/Input-File-Format  #
# for a thorough description of the format and its components.             #
#                                                                          #
############################################################################

# Declare special types.
@enum FLOW_DIRECTION POSITIVE=1 NEGATIVE=-1 UNKNOWN=0

function parse_epanet_file(path::String)
    file_contents = readstring(open(path))
    file_contents = replace(file_contents, "\t", "    ")
    lines = split(file_contents, '\n')

    section = headings = nothing
    headings_exist = false
    epanet_dict = Dict{String, Any}()

    for (i, line) in enumerate(lines)
        if ismatch(r"^\s*\[(.*)\]", line) # If a section heading.
            section = lowercase(strip(line, ['[', ']']))
            headings_exist = ismatch(r"^;", lines[i+1])
            if !headings_exist
                epanet_dict["$section"] = Dict{String, Any}()
            end
        elseif ismatch(r"^;", line) # If a section heading.
            headings = split(lowercase(strip(line, [';'])))
            epanet_dict["$section"] = Dict{String, Array}(h => [] for h = headings)
        elseif length(line) == 0
            continue
        else
            if headings_exist
                data = split(lowercase(strip(line, [';'])))
                for (j, heading) in enumerate(headings)
            if j <= length(data)
                push!(epanet_dict["$section"]["$heading"], data[j])
            else
                push!(epanet_dict["$section"]["$heading"], "")
            end
                end
            else
                heading = strip(split(line, "  ")[1])
                data = split(lowercase(strip(replace(line, heading, ""), [';'])))
                epanet_dict["$section"][lowercase(heading)] = data
            end
        end
    end

    # Parse relevant data into a more structured format.
    dict = Dict{String, Any}()
    dict["title"] = parse_title(epanet_dict["title"])
    dict["options"] = parse_options(epanet_dict["options"])
    dict["junctions"] = parse_junctions(epanet_dict["junctions"], dict["options"])
    dict["pipes"] = parse_pipes(epanet_dict["pipes"], dict["options"])
    dict["reservoirs"] = parse_reservoirs(epanet_dict["reservoirs"], dict["options"])
    dict["tanks"] = parse_tanks(epanet_dict["tanks"])
    dict["valves"] = parse_valves(epanet_dict["valves"])
    dict["multinetwork"] = false
    dict["per_unit"] = false

    return dict
end

function allequal(x) 
    return all(y->y == x[1], x)
end

function parse_general(dtype::Type, data::Any)
    do_not_parse = dtype == String || dtype == FLOW_DIRECTION || dtype == typeof(data)
    return do_not_parse ? data : parse(dtype, data)
end

function parse_title(data::Dict{String, Any})
    return length(keys(data)) > 0 ? first(keys(data)) : ""
end

function parse_junctions(data::Dict{String, Array}, options::Dict{String, Any})
    # Get the demand units (e.g., LPS, GPM).
    demand_units = options["units"]

    # Initialize scalars to convert data to SI units.
    demand_scalar = nothing
    elev_scalar = nothing

    if demand_units == "lps" # If liters per second...
        # Convert from liters per second to cubic meters per second.
        demand_scalar = 1.0e-3 * options["demand_multiplier"]

        # Retain the original value (in meters).
        elev_scalar = 1.0
    elseif demand_units == "gpm" # If gallons per minute...
        # Convert from gallons per minute to cubic meters per second.
        demand_scalar = 6.30902e-5 * options["demand_multiplier"]

        # Convert elevation from feet to meters.
        elev_scalar = 0.3048
    else
        error("Could not find a valid \"units\" option type.")
    end

    # Convert demand data to SI units.
    data["demand"] = [parse(Float64, demand) for demand in data["demand"]]
    data["demand"] = demand_scalar .* data["demand"]

    # Convert elev data to SI units.
    data["elev"] = [parse(Float64, elev) for elev in data["elev"]]
    data["elev"] = elev_scalar .* data["elev"]

    # Specify the data types for the junction data.
    columns = Dict("demand" => Float64, "elev" => Float64,
                   "id" => String, "pattern" => String)

    # Ensure the arrays describing junction data are all of equal lengths.
    @assert(allequal([length(data[column]) for column in keys(columns)]))

    # Return an array of junction dictionaries with the correct data types.
    arr = [Dict(c => parse_general(v, data[c][i]) for (c, v) in columns) for i = 1:length(data["id"])]

    # Scale the quantities appropriately.
    return Dict{String, Any}(data["id"][i] => arr[i] for i = 1:length(arr))
end

function parse_pipes(data::Dict{String, Array}, options::Dict{String, Any})
    # Get the demand units (e.g., LPS, GPM).
    demand_units = options["units"]

    # Get the headloss type (i.e., Darcy-Weisbach or Hazen-Williams)
    headloss_type = options["headloss"]

    # Initialize scalars to convert data to SI units.
    diameter_scalar = nothing
    length_scalar = nothing
    roughness_scalar = nothing

    if demand_units == "lps" # If liters per second...
        # Convert diameter from millimeters to meters.
        diameter_scalar = 0.001

        # Retain the original value (in meters).
        length_scalar = 1.0

        if headloss_type == "d-w"
            # Convert roughness from millimeters to meters.
            roughness_scalar = 0.001
        elseif headloss_type == "h-w"
            # Retain the original value (unitless).
            roughness_scalar = 1.0
        end
    elseif demand_units == "gpm" # If gallons per minute...
        # Convert diameter from inches to meters.
        diameter_scalar = 0.0254

        # Convert length from feet to meters.
        length_scalar = 0.3048

        if headloss_type == "d-w"
            # Convert roughness from millifeet to meters.
            roughness_scalar = 3.048e-4
        elseif headloss_type == "h-w"
            # Retain the original value (unitless).
            roughness_scalar = 1.0
        end
    else
        error("Could not find a valid \"units\" option type.")
    end

    # Convert diameter data to SI units.
    data["diameter"] = [parse(Float64, diameter) for diameter in data["diameter"]]
    data["diameter"] = diameter_scalar .* data["diameter"]

    # Convert length data to SI units.
    data["length"] = [parse(Float64, length) for length in data["length"]]
    data["length"] = length_scalar .* data["length"]

    # Convert roughness data to SI units.
    data["roughness"] = [parse(Float64, roughness) for roughness in data["roughness"]]
    data["roughness"] = roughness_scalar .* data["roughness"]

    # Specify the data types for the pipe data.
    columns = Dict("diameter" => Float64, "id" => String, "length" => Float64,
                   "minorloss" => Float64, "node1" => String, "node2" => String,
                   "roughness" => Float64, "status" => String,
                   "flow_direction" => FLOW_DIRECTION)

    # Populate the flow direction data.
    data["flow_direction"] = Array{FLOW_DIRECTION}(length(data["id"]))
    fill!(data["flow_direction"], UNKNOWN) # The initial flow direction is unknown.

    # Ensure the arrays describing pipe data are all of equal lengths.
    @assert(allequal([length(data[column]) for column in keys(columns)]))

    # Create a dictionary of pipe dictionaries with the correct data types.
    arr = [Dict(c => parse_general(v, data[c][i]) for (c, v) in columns) for i = 1:length(data["id"])]
    return Dict{String, Any}(data["id"][i] => arr[i] for i = 1:length(arr))
end

function parse_reservoirs(data::Dict{String, Array}, options::Dict{String, Any})
    demand_units = options["units"]
    head_scalar = nothing

    if demand_units == "lps" # If liters per second...
        # Retain the original value (in meters).
        head_scalar = 1.0
    elseif demand_units == "gpm" # If gallons per minute...
        # Convert from feet to meters.
        head_scalar = 0.3048
    else
        error("Could not find a valid \"units\" option type.")
    end

    # Convert diameter data to SI units.
    data["head"] = [parse(Float64, head) for head in data["head"]]
    data["head"] = head_scalar .* data["head"]

    # Specify the data types for the reservoir data.
    columns = Dict("head" => Float64, "id" => String, "pattern" => String)

    # Ensure the arrays describing reservoir data are all of equal lengths.
    @assert(allequal([length(data[column]) for column in keys(columns)]))

    # Return an array of reservoir dictionaries with the correct data types.
    arr = [Dict(c => parse_general(v, data[c][i]) for (c, v) in columns) for i = 1:length(data["id"])]
    return Dict{String, Any}(data["id"][i] => arr[i] for i = 1:length(arr))
end

function parse_tanks(data::Dict{String, Array})
    # Specify the data types for the reservoir data.
    columns = Dict("diameter" => Float64, "elevation" => Float64,
                   "id" => String, "initlevel" => Float64,
                   "maxlevel" => Float64, "minlevel" => Float64,
                   "minvol" => Float64, "volcurve" => String)

    # Ensure the arrays describing reservoir data are all of equal lengths.
    @assert(allequal([length(data[column]) for column in keys(columns)]))

    # Return an array of reservoir dictionaries with the correct data types.
    arr = [Dict(c => parse_general(v, data[c][i]) for (c, v) in columns) for i = 1:length(data["id"])]
    return Dict{String, Any}(data["id"][i] => arr[i] for i = 1:length(arr))
end

function parse_valves(data::Dict{String, Array})
    # Specify the data types for the reservoir data.
    columns = Dict("diameter" => Float64, "id" => String,
                   "minorloss" => Float64, "node1" => String,
                   "node2" => String, "setting" => Float64,
                   "type" => String)

    # Ensure the arrays describing reservoir data are all of equal lengths.
    @assert(allequal([length(data[column]) for column in keys(columns)]))

    # Return an array of reservoir dictionaries with the correct data types.
    arr = [Dict(c => parse_general(v, data[c][i]) for (c, v) in columns) for i = 1:length(data["id"])]
    return Dict{String, Any}(data["id"][i] => arr[i] for i = 1:length(arr))
end

function parse_options(data::Dict{String, Any})
    units = data["units"][1]
    headloss = data["headloss"][1]
    demand_multiplier = parse(Float64, data["demand multiplier"][1])
    viscosity = parse(Float64, data["viscosity"][1]) * 1.0e-3
    return Dict{String, Any}("units" => units, "headloss" => headloss,
                             "demand_multiplier" => demand_multiplier,
                             "viscosity" => viscosity)
end
