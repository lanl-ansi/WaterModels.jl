function _initialize_epanet_dictionary()
    data = Dict{String,Any}("time_series" => Dict{String,Any}(), "top_comments" => [])
    data["section"] = Dict{String,Array}(x => [] for x in _INP_SECTIONS)
    return data # Return an empty dictionary of the sectional data.
end


function _read_epanet_sections(file_path::String)
    # Initialize the dictionary that stores raw EPANET data.
    data = _initialize_epanet_dictionary()

    # Instantiate metadata used in parsing EPANET sections.
    section, line_number = nothing, 0

    # Populate sections in the EPANET dictionary.
    for line in eachline(file_path)
        # Update important data associated with the line.
        line, line_number = strip(line), line_number + 1

        if length(line) == 0
            # If the line is empty, continue to the next line.
            continue
        elseif startswith(line, "[")
            # If the line starts with "[", it must be a section.
            values = split(line; limit = 1)
            section_tmp = uppercase(values[1])

            if section_tmp in _INP_SECTIONS
                # If the section title is valid, store it as the current.
                section = section_tmp
                continue
            elseif section_tmp == "[END]"
                # If the section title is "[END]", parsing is complete.
                break
            else
                # If the section title is not valid, throw an error.
                msg = "File \"$(file_path)\", line $(line_number) has invalid section \"$(section_tmp)\"."
                Memento.error(_LOGGER, msg)
            end
        elseif section === nothing && startswith(line, ";")
            # If there is no current section and the section is a comment, store it.
            data["top_comments"] = vcat(data["top_comments"], line[2:end])
            continue
        elseif section === nothing && !startswith(line, ";")
            # If there is no current section and the line is not a comment, throw an error.
            msg = "File \"$(file_path)\", line $(line_number) has invalid syntax."
            Memento.error(_LOGGER, msg)
        end

        # Append a tuple of the line number and the line string to the dictionary.
        push!(data["section"][section], (line_number, line))
    end

    return data # Return the "raw" EPANET data dictionary.
end


"""
    parse_epanet(path)

Parses an [EPANET](https://www.epa.gov/water-research/epanet) (.inp) file from
the file path `path` and returns a WaterModels data structure (a dictionary of
data). See the [OpenWaterAnalytics
Wiki](https://github.com/OpenWaterAnalytics/EPANET/wiki/Input-File-Format) for
a thorough description of the EPANET format and its components.
"""
function parse_epanet(filename::String)
    # Read in raw sectional EPANET data.
    epanet_data = _read_epanet_sections(filename)

    # Parse [OPTIONS] section.
    _parse_epanet_options(epanet_data)

    # Parse [TIMES] section.
    _read_time!(epanet_data)

    # Parse [CURVES] section.
    _read_curve!(epanet_data)

    # Parse [PATTERNS] section.
    _read_pattern!(epanet_data)

    # Parse [JUNCTIONS] section.
    _read_demand!(epanet_data)

    # Parse [RESERVOIRS] section.
    _read_reservoir!(epanet_data)

    # Parse [TANKS] section.
    _read_tank!(epanet_data)

    # Add nodes for components with nodal properties.
    _add_nodes!(epanet_data)

    # Add a data structure mapping names to node indices.
    _add_node_map!(epanet_data)

    # Update time series information for reservoirs.
    _update_reservoir_ts!(epanet_data)

    # Parse [PIPES] section.
    _read_pipe!(epanet_data)

    # Set up a dictionary containing valve objects.
    epanet_data["valve"] = Dict{String,Any}()

    # Set up a dictionary containing short pipe objects.
    epanet_data["short_pipe"] = Dict{String,Any}()

    # Set up a dictionary containing expansion short pipe objects.
    epanet_data["ne_short_pipe"] = Dict{String,Any}()

    # Set up a dictionary containing expansion pump objects.
    epanet_data["ne_pump"] = Dict{String,Any}()

    # Set up a dictionary containing design pipe objects.
    epanet_data["des_pipe"] = Dict{String,Any}()

    # Parse [PUMPS] section.
    _read_pump!(epanet_data)

    # Parse [VALVES] section.
    _read_regulator!(epanet_data)

    # Parse [ENERGY] section.
    _read_energy!(epanet_data)

    # Parse [STATUS] section.
    _read_status!(epanet_data)

    # Parse [CONTROLS] section.
    _read_controls!(epanet_data)

    # Update pump data using data from the [ENERGY] section.
    _update_pump_energy!(epanet_data)

    # Parse [COORDINATES] section.
    _read_coordinate!(epanet_data)

    # Parse [TITLE] section.
    _read_title!(epanet_data)

    # Add or remove time series information.
    _update_time_series!(epanet_data)

    # Delete redundant data that may have been parsed or added to the dictionary.
    for key in [
        "section",
        "node_map",
        "curve",
        "pattern",
        "hydraulic_time_step",
        "pattern_time_step",
        "demand_multiplier",
        "flow_units",
        "duration",
        "top_comments",
        "default_pattern",
        "energy_price",
        "energy_pattern",
        "energy_efficiency",
        "duration",
        "start_time",
        "demand_charge",
    ]
        key in keys(epanet_data) && delete!(epanet_data, key)
    end

    # Add other data required by InfrastructureModels.
    epanet_data["per_unit"] = false

    # Remove all keys that have values of nothing.
    _clean_nothing!(epanet_data)

    # Return the dictionary.
    return epanet_data
end


function _update_time_series!(data::Dict{String,<:Any})
    duration, time_step = data["duration"], data["hydraulic_time_step"]
    num_steps = convert(Int, floor(duration * inv(time_step))) + 1
    data["time_step"] = convert(Float64, time_step)

    if num_steps >= 1 && keys(data["pattern"]) != ["1"]
        for type in ["demand", "node", "pump"]
            length(data["time_series"][type]) == 0 && delete!(data["time_series"], type)
        end

        data["time_series"]["duration"] = data["duration"]
        data["time_series"]["time_step"] = data["hydraulic_time_step"]
        data["time_series"]["num_steps"] = num_steps

        for pattern in keys(data["pattern"])
            if !(length(data["pattern"][pattern]) in [1, num_steps - 1])
                Memento.error(_LOGGER, "Pattern \"$(pattern)\" is of the wrong length.")
            end
        end
    else
        delete!(data, "time_series")
    end
end

"Standardize component indices to be integers instead of strings."
function _transform_component_indices(components::Dict{String,<:Any})
    if all([tryparse(Int, x["name"]) !== nothing for (i, x) in components])
        for (name, component) in components
            component["index"] = parse(Int, component["name"])
        end
    else
        for (name, component) in components
            component["index"] = parse(Int, component["index"])
        end
    end

    return Dict{String,Any}(string(x["index"]) => x for (i, x) in components)
end


"Standardize component time series indices to match primary indices."
function _transform_time_series_indices(components::Dict{String,<:Any}, time_series::Dict{String,<:Any})
    name_to_index = Dict{String,String}(x["name"] => i for (i, x) in components)
    return Dict{String,Any}(name_to_index[i] => x for (i, x) in time_series)
end


function _update_demand_ts!(data::Dict{String,<:Any})
    # Create a temporary dictionary representing the demand time series.
    demand_ts = Dict{String,Any}()

    # Ensure that demand time series data use new demand indices.
    for (i, demand) in data["time_series"]["demand"]
        key = findfirst(x -> x["source_id"][2] == i, data["demand"])
        demand_ts[key] = demand
    end

    # Update the demand entry in the time series dictionary.
    return data["time_series"]["demand"] = demand_ts
end

function _update_reservoir_ts!(data::Dict{String,<:Any})
    # Create a temporary dictionary representing the reservoir time series.
    node_ts = Dict{String,Any}()

    # Ensure that reservoir time series data use new reservoir indices.
    for (i, reservoir) in data["time_series"]["reservoir"]
        key = findfirst(x -> x["source_id"][2] == i, data["reservoir"])
        node_index = data["reservoir"][key]["node"]
        node_ts[string(node_index)] = reservoir
    end

    # Update the reservoir entry in the time series dictionary.
    data["time_series"]["node"] = node_ts
    return delete!(data["time_series"], "reservoir")
end

function _add_nodes!(data::Dict{String,<:Any})
    # Define EPANET nodal types that should be attached to nodes.
    comp_types = ["demand", "reservoir", "tank"]

    # Obtain the original EPANET indices for all nodal components.
    comp_names = vcat([collect(keys(data[t])) for t in comp_types]...)

    # Set a flag based on if there are duplicates in comp_names.
    keep_name = length(unique(comp_names)) == length(comp_names)

    # Initialize the node dictionary and temporary index.
    data["node"], index = Dict{String,Any}(), 0

    # Loop over EPANET nodal types.
    for comp_type in comp_types
        # Loop over all nodal components of a particular type.
        for (comp_name, comp) in data[comp_type]
            # Initialize the node entry to be added.
            key = keep_name ? comp_name : string(index += 1)
            node = Dict{String,Any}("name" => key, "index" => parse(Int, key))

            # Save common data from the component within the node entry.
            node["status"], node["elevation"] = comp["status"], pop!(comp, "elevation")

            if haskey(node, "head_nominal")
                node["head_nominal"] = pop!(comp, "head_nominal")
            else node["elevation"]
                node["head_nominal"] = node["elevation"]
            end

            node["source_id"] = comp["source_id"]

            # Append the node to the node dictionary.
            data["node"][key] = node

            # Store the index of the node in the component dictionary.
            comp["node"] = node["index"]
        end
    end
end


function _split_line(line::AbstractString)
    _vc = split(line, ","; limit = 1)
    _values = nothing
    _comment = nothing

    if length(_vc) != 0
        if length(_vc) == 1
            _values = split(_vc[1])
        elseif _vc[1] == ""
            _comment = _vc[2]
        else
            _values = split(_vc[1])
            _comment = _vc[2]
        end
    end

    return _values, _comment
end


"""
Converts EPANET clocktime format to seconds.
Parameters
----------
s : string
    EPANET time string. Options are 'HH:MM:SS', 'HH:MM', HH'
am : string
    options are AM or PM
Returns
-------
Integer value of time in seconds
"""
function _clock_time_to_seconds(s::AbstractString, am_pm::AbstractString)
    if uppercase(am_pm) == "AM"
        am = true
    elseif uppercase(am_pm) == "PM"
        am = false
    else
        Memento.error(_LOG, "AM/PM option not recognized (options are AM or PM")
    end

    time_tuple = match(r"^(\d+):(\d+):(\d+)$", s)

    if time_tuple !== nothing
        seconds_1 = parse(Int, time_tuple[1]) * 3600
        seconds_2 = parse(Int, time_tuple[2]) * 60
        seconds_3 = convert(Int, round(parse(Float64, time_tuple[3])))
        time_seconds = seconds_1 + seconds_2 + seconds_3

        if startswith(s, "12")
            time_seconds -= 3600 * 12
        end

        if !am
            if time_seconds >= 3600 * 12
                Memento.error(_LOG, "Cannot specify AM/PM for times greater than 12:00:00")
            end

            time_seconds += 3600 * 12
        end

        return time_seconds
    else
        time_tuple = match(r"^(\d+):(\d+)$", s)

        if time_tuple !== nothing
            seconds_1 = parse(Int, time_tuple[1]) * 3600
            seconds_2 = parse(Int, time_tuple[2]) * 60

            if startswith(s, "12")
                time_seconds -= 3600 * 12
            end

            if !am
                if time_seconds >= 3600 * 12
                    Memento.error(
                        _LOG,
                        "Cannot specify AM/PM for times greater than 12:00:00",
                    )
                end

                time_seconds += 3600 * 12
            end

            return time_seconds
        else
            time_tuple = match(r"^(\d+)$", s)

            if time_tuple !== nothing
                time_seconds = parse(Int, time_tuple[1]) * 3600

                if startswith(s, "12")
                    time_seconds -= 3600 * 12
                end

                if !am
                    if time_seconds >= 3600 * 12
                        Memento.error(
                            _LOG,
                            "Cannot specify AM/PM for times greater than 12:00:00",
                        )
                    end

                    time_seconds += 3600 * 12
                end

                return time_seconds
            else
                Memento.error(_LOGGER, "Time format in INP file not recognized")
            end
        end
    end
end

"""
Converts EPANET time format to seconds.
Parameters
----------
s : string
    EPANET time string. Options are 'HH:MM:SS', 'HH:MM', 'HH'
Returns
-------
 Integer value of time in seconds.
"""
function _string_time_to_seconds(s::AbstractString)
    time_tuple = match(r"^(\d+):(\d+):(\d+)$", s)

    if time_tuple !== nothing
        seconds_1 = parse(Int, time_tuple[1]) * 3600
        seconds_2 = parse(Int, time_tuple[2]) * 60
        seconds_3 = convert(Int, round(parse(Float64, time_tuple[3])))
        return seconds_1 + seconds_2 + seconds_3
    else
        time_tuple = match(r"^(\d+):(\d+)$", s)

        if time_tuple !== nothing
            seconds_1 = parse(Int, time_tuple[1]) * 3600
            seconds_2 = parse(Int, time_tuple[2]) * 60
            return seconds_1 + seconds_2
        else
            time_tuple = match(r"^(\d+)$", s)

            if time_tuple !== nothing
                return parse(Int, time_tuple[1]) * 3600
            else
                Memento.error(_LOGGER, "Time format in INP file not recognized")
            end
        end
    end
end

function _clean_nothing!(data)
    for (key, value) in data
        if isa(value, Dict)
            value = _clean_nothing!(value)
        elseif value === nothing
            delete!(data, key)
        end
    end
end

function _read_coordinate!(data::Dict{String,<:Any})
    for (line_number, line) in data["section"]["[COORDINATES]"]
        current = split(split(line, ";")[1])
        length(current) == 0 && continue # Skip.

        # Add coordinate data to the corresponding node.
        node_id = data["node_map"][current[1]]
        x, y = parse.(Float64, current[2:3])
        data["node"][string(node_id)]["coordinates"] = (x, y)
    end
end

function _read_curve!(data::Dict{String,<:Any})
    # Initialize dictionary associated with curves.
    data["curve"] = Dict{String,Array{Tuple{Float64,Float64}}}()

    # Loop over all lines in the [CURVES] section and parse each.
    for (line_number, line) in data["section"]["[CURVES]"]
        current = split(split(line, ";")[1])
        length(current) == 0 && continue # Skip.
        curve_name = current[1]

        # If `curve_name` is not in the dictionary, initialize it.
        if !(curve_name in keys(data["curve"]))
            data["curve"][curve_name] = Array{Tuple{Float64,Float64},1}()
        end

        # Add the point specified on the line the curve entry.
        x, y = parse(Float64, current[2]), parse(Float64, current[3])
        data["curve"][curve_name] = vcat(data["curve"][curve_name], (x, y))
    end
end

function _read_energy!(data::Dict{String,<:Any})
    # Loop over all lines in the [ENERGY] section and parse each.
    for (line_number, line) in data["section"]["[ENERGY]"]
        current = split(split(line, ";")[1])
        length(current) == 0 && continue # Skip.

        if uppercase(current[1]) == "GLOBAL"
            # Read global energy data (e.g., price, default pump efficiency).
            if uppercase(current[2]) == "PRICE"
                data["energy_price"] = inv(3.6e6) * parse(Float64, current[3])
            elseif uppercase(current[2]) == "PATTERN"
                data["energy_pattern"] = current[3]
            elseif uppercase(current[2]) in ["EFFIC", "EFFICIENCY"]
                data["energy_efficiency"] = 0.01 * parse(Float64, current[3])
            else
                Memento.warn(_LOGGER, "Unknown entry in ENERGY section: $(line)")
            end
        elseif uppercase(current[1]) == "DEMAND"
            # Read in the demand charge, an additional cost per maximum usage.
            data["demand_charge"] = parse(Float64, current[3])
        elseif uppercase(current[1]) == "PUMP"
            # Get the corresponding pump data object.
            pump_key = findfirst(x -> current[2] == x["source_id"][2], data["pump"])
            pump = data["pump"][pump_key]

            if uppercase(current[3]) == "PRICE"
                # Parse the cost of pump operation.
                price = parse(Float64, current[4]) # Price per kilowatt hour.
                pump["energy_price"] = price * inv(3.6e6) # Price per Joule.
            elseif uppercase(current[3]) == "PATTERN"
                # Read in the pattern for scaling the pump's energy price.
                pump["energy_pattern"] = current[4]
            elseif uppercase(current[3]) in ["EFFIC", "EFFICIENCY"]
                # Obtain and scale head-versus-flow pump curve.
                x = first.(data["curve"][current[4]]) # Flow rate.
                y = 0.01 .* last.(data["curve"][current[4]]) # Efficiency.

                if data["flow_units"] == "LPS" # If liters per second...
                    # Convert from liters per second to cubic meters per second.
                    x *= 1.0e-3
                elseif data["flow_units"] == "CMH" # If cubic meters per hour...
                    # Convert from cubic meters per hour to cubic meters per second.
                    x *= inv(3600.0)
                elseif data["flow_units"] == "GPM" # If gallons per minute...
                    # Convert from gallons per minute to cubic meters per second.
                    x *= 6.30902e-5
                else
                    Memento.error(_LOGGER, "Could not find a valid \"units\" option type.")
                end

                for k = 1:length(y)
                    y[k] = x[k] > 0.0 && y[k] > 1.0e-2 ? y[k] : 1.0e-2
                end

                # Curve of efficiency (unitless) versus the flow rate through a pump.
                pump["efficiency_curve"] = Array([(x[j], y[j]) for j = 1:length(x)])
            else
                Memento.warn(_LOGGER, "Unknown entry in ENERGY section: $(line)")
            end
        else
            Memento.warn(_LOGGER, "Unknown entry in ENERGY section: $(line)")
        end
    end
end

function _read_status!(data::Dict{String,<:Any})
    # Loop over all lines in the [STATUS] section and parse each.
    for (line_number, line) in data["section"]["[STATUS]"]
        current = split(split(line, ";")[1])
        length(current) == 0 && continue # Skip.
        component_id = current[1]

        if haskey(data["pipe"], component_id)
            # Ensure the pipe is categorized as having a valve.
            data["pipe"][component_id]["has_valve"] = true
        end
    end
end

function _read_controls!(data::Dict{String,<:Any})
    # Loop over all lines in the [CONTROLS] section and parse each.
    for (line_number, line) in data["section"]["[CONTROLS]"]
        current = split(split(line, ";")[1])
        length(current) == 0 && continue # Skip.
        component_id = current[2]

        if uppercase(current[1]) == "LINK" && haskey(data["pipe"], component_id)
            data["pipe"][component_id]["has_valve"] = true
        end
    end
end

function _update_pump_energy!(data::Dict{String,<:Any})
    # Use data previously parsed to infer missing pump data.
    for (pump_id, pump) in data["pump"]
        # If the pump does not have an efficiency curve, assume a default value.
        if !("efficiency_curve" in keys(pump)) && "energy_efficiency" in keys(data)
            # If the pump does not have an efficiency, assume the global option.
            pump["efficiency"] = data["energy_efficiency"]
        elseif !("efficiency_curve" in keys(pump))
            # Otherwise, assume the pump is perfectly efficient.
            Memento.warn(
                _LOGGER,
                "Efficiency for pump \"$(pump["name"])\" could not be found.",
            )
            pump["efficiency"] = 1.0
        end

        # If the pump does not have an energy pattern, assume the global option.
        if !("energy_pattern" in keys(pump)) && "energy_pattern" in keys(data)
            pump["energy_pattern"] = data["energy_pattern"]
        end

        # If the energy price is not specified for the pump, set a reasonable value.
        if !("energy_price" in keys(pump)) && "energy_price" in keys(data)
            # If the pump does not have an energy price, assume the global option.
            pump["energy_price"] = data["energy_price"]
        elseif !("energy_price" in keys(pump)) # If there is no energy price...
            # Otherwise, assume the pump does not have an energy price.
            Memento.warn(_LOGGER, "Price for pump \"$(pump["name"])\" could not be found.")
            pump["energy_price"] = 0.0
        end

        # If an energy pattern is specified, store the time series of prices.
        if "energy_pattern" in keys(pump)
            # Build the pattern of prices using existing data.
            pattern_name = pump["energy_pattern"]
            price_pattern = pump["energy_price"] .* data["pattern"][pattern_name]
            push!(price_pattern, minimum(price_pattern))

            # Add the varying price data to the pump time series entry.
            entry = Dict{String,Array{Float64}}("energy_price" => price_pattern)
            data["time_series"]["pump"][pump_id] = entry

            # Delete the "energy_pattern" data from the pump.
            delete!(pump, "energy_pattern")
        end
    end
end

"Parse and store demand data from an EPANET-formatted data dictionary."
function _read_demand!(data::Dict{String,<:Any})
    # Initialize dictionaries associated with demands.
    data["demand"] = Dict{String,Dict{String,Any}}()
    data["time_series"]["demand"] = Dict{String,Any}()

    # Initialize a temporary index to be updated while parsing.
    index::Int = 0 # This will track the order that demands appear.

    # Loop over all lines in the [JUNCTIONS] section and parse each.
    for (line_number, line) in data["section"]["[JUNCTIONS]"]
        current = split(split(line, ";")[1])
        length(current) == 0 && continue # Skip.

        # Initialize the demand entry to be added.
        demand = Dict{String,Any}("name" => current[1], "status" => STATUS_ACTIVE)
        demand["source_id"] = ["demand", current[1]]
        demand["dispatchable"] = false

        # Parse the elevation of the demand (in meters).
        if data["flow_units"] == "LPS" || data["flow_units"] == "CMH"
            # Retain the original value (in meters).
            demand["elevation"] = parse(Float64, current[2])
        elseif data["flow_units"] == "GPM" # If gallons per minute...
            # Convert elevation from feet to meters.
            demand["elevation"] = 0.3048 * parse(Float64, current[2])
        else
            Memento.error(_LOGGER, "Could not find a valid \"units\" option type.")
        end

        # Parse the demand at the demand (in cubic meters per second).
        if length(current) > 2 # A value for demand must exist.
            # Calculate the unscaled (unit unknown) demand.
            demand_unscaled = parse(Float64, current[3]) * data["demand_multiplier"]

            # Convert the unscaled demand to cubic meters per second.
            if data["flow_units"] == "LPS" # If liters per second...
                # Convert from liters per second to cubic meters per second.
                demand["flow_nominal"] = 1.0e-3 * demand_unscaled
            elseif data["flow_units"] == "CMH" # If cubic meters per hour...
                # Convert from cubic meters per hour to cubic meters per second.
                demand["flow_nominal"] = inv(3600.0) * demand_unscaled
            elseif data["flow_units"] == "GPM" # If gallons per minute...
                # Convert from gallons per minute to cubic meters per second.
                demand["flow_nominal"] = 6.30902e-5 * demand_unscaled
            end
        else # No value for demand exists in the line.
            demand["status"] = STATUS_INACTIVE # Set the status of the demand to zero.
            demand["flow_nominal"] = 0.0 # Initialize demand to zero.
        end

        demand["flow_min"] = demand["flow_nominal"]
        demand["flow_max"] = demand["flow_nominal"]

        # Parse the name of the pattern used to scale demand at the demand.
        pattern = length(current) > 3 ? current[4] : data["default_pattern"]

        # Scale the parsed demand by the specified pattern, if it exists.
        if pattern !== nothing && length(data["pattern"][pattern]) > 1
            flow_rate = demand["flow_nominal"] .* data["pattern"][pattern]
            push!(flow_rate, minimum(flow_rate))

            entry = Dict{String,Array{Float64}}("flow_nominal" => flow_rate)
            entry["flow_min"] = entry["flow_nominal"]
            entry["flow_max"] = entry["flow_nominal"]
            data["time_series"]["demand"][current[1]] = entry

            demand["flow_min"] = minimum(entry["flow_min"])
            demand["flow_max"] = maximum(entry["flow_max"])
        elseif pattern !== nothing && pattern != "1"
            demand["flow_nominal"] *= data["pattern"][pattern][1]
            demand["flow_min"] = demand["flow_nominal"]
            demand["flow_max"] = demand["flow_nominal"]
        end

        # Add a temporary index to be used in the data dictionary.
        demand["index"] = string(index += 1)

        # Append the demand entry to the data dictionary.
        data["demand"][current[1]] = demand
    end

    # Replace with new dictionaries that use integer component indices.
    data["demand"] = _transform_component_indices(data["demand"])
    data["time_series"]["demand"] = _transform_time_series_indices(
        data["demand"], data["time_series"]["demand"])
end

function _add_node_map!(data::Dict{String,<:Any})
    # Initialize the dictionary mapping component names to node indices.
    data["node_map"] = Dict{String,Int}()

    # Map component names to node indices.
    for type in ["demand", "reservoir", "tank"]
        map = Dict{String,Int}(x["name"] => x["node"] for (i, x) in data[type])
        data["node_map"] = merge(data["node_map"], map)
    end

    # Ensure the number of node-type components is equal to the number of nodes.
    num_components = sum(length(data[type]) for type in ["demand", "reservoir", "tank"])

    # Ensure the number of node-type components is equal to the length of the map.
    if num_components != length(data["node_map"])
        Memento.error(_LOGGER, "Duplicate node-type component name detected.")
    end
end

function _parse_epanet_options(data::Dict{String,<:Any})
    # Loop over all lines in the [OPTIONS] section and parse each.
    for (line_number, line) in data["section"]["[OPTIONS]"]
        words, comments = _split_line(line)

        if words !== nothing && length(words) > 0
            if length(words) < 2
                Memento.error(_LOGGER, "No value provided for option $(words[1])")
            end

            key = uppercase(words[1])

            if key == "UNITS"
                data["flow_units"] = uppercase(words[2])
            elseif key == "HEADLOSS"
                data["head_loss"] = uppercase(words[2])
            elseif key == "VISCOSITY"
                data["viscosity"] = 1.0e-3 * parse(Float64, words[2])
            elseif key == "PATTERN"
                data["default_pattern"] = words[2]
            elseif key == "DEMAND"
                if length(words) > 2
                    data["demand_multiplier"] = parse(Float64, words[3])
                else
                    Memento.error(_LOGGER, "No value provided for Demand Multiplier")
                end
            end
        end
    end

    # Ensure that required data fields have been populated.
    @assert all(haskey.(Ref(data), ["flow_units", "head_loss", "viscosity"]))
end

function _read_pattern!(data::Dict{String,<:Any})
    # Initialize dictionary associated with patterns.
    data["pattern"] = Dict{String,Array}()

    # Loop over all lines in the [PATTERNS] section and parse each.
    for (line_number, line) in data["section"]["[PATTERNS]"]
        current = split(split(line, ";")[1])
        length(current) == 0 && continue # Skip.
        name = current[1]

        if !(name in keys(data["pattern"]))
            data["pattern"][name] = Vector{Float64}()
        end

        values = [parse(Float64, s) for s in current[2:end]]
        data["pattern"][name] = vcat(data["pattern"][name], values)
    end

    # Ensure the hydraulic and pattern time steps are equal.
    if "pattern_time_step" in keys(data)
        # Get both time steps defined in the data.
        hydraulic_time_step = data["hydraulic_time_step"]
        pattern_time_step = data["pattern_time_step"]

        if hydraulic_time_step != pattern_time_step
            # Compute how many pattern steps are in one hydraulic step.
            factor = convert(Int, pattern_time_step / hydraulic_time_step)

            # For all patterns, ensure they match the hydraulic periods.
            for (pattern_name, pattern) in data["pattern"]
                pattern = [pattern[div(i, factor)+1] for i = 0:(factor*length(pattern)-1)]
                data["pattern"][pattern_name] = pattern
            end
        end
    end

    # Get the number of steps that should now exist in each pattern.
    duration, time_step = data["duration"], data["hydraulic_time_step"]
    num_steps = convert(Int, floor(duration / time_step))

    # Ensure patterns of length one are equal to the real pattern length.
    for (pattern_name, pattern) in data["pattern"]
        if length(pattern) == 1 && num_steps != 1
            data["pattern"][pattern_name] = ones(num_steps) * pattern[1]
        end
    end

    if data["default_pattern"] == "" && "1" in keys(data["pattern"])
        # If there is a pattern called "1", then it is the default pattern if no other is supplied.
        data["default_pattern"] = "1"
    elseif !(data["default_pattern"] in keys(data["pattern"]))
        # Sanity check: if the default pattern does not exist and it is not "1", then balk.
        # If default is "1" but it does not exist, then it is constant.
        # Any other default that does not exist is an error.
        if data["default_pattern"] !== nothing && data["default_pattern"] != "1"
            Memento.error("Default pattern \"$(data["default_pattern"])\" is undefined.")
        end

        data["default_pattern"] = nothing
    end
end

function _read_pipe!(data::Dict{String,<:Any})
    # Initialize dictionary associated with pipes.
    data["pipe"] = Dict{String,Dict{String,Any}}()

    # Initialize a temporary index to be updated while parsing.
    index::Int = 0

    # Loop over all lines in the [PIPES] section and parse each.
    for (line_number, line) in data["section"]["[PIPES]"]
        current = split(split(line, ";")[1])
        length(current) == 0 && continue # Skip.

        # Initialize the pipe entry to be added.
        pipe = Dict{String,Any}("name" => current[1], "status" => STATUS_ACTIVE)
        pipe["source_id"] = ["pipe", current[1]]
        pipe["node_fr"] = data["node_map"][current[2]]
        pipe["node_to"] = data["node_map"][current[3]]

        # Store all measurements associated with pipes in metric units.
        if data["flow_units"] == "LPS" || data["flow_units"] == "CMH"
            # Retain the original value (in meters).
            pipe["length"] = parse(Float64, current[4])

            # Convert diameter from millimeters to meters.
            pipe["diameter"] = 0.001 * parse(Float64, current[5])

            if data["head_loss"] == "D-W" # If Darcy-Weisbach head loss is used...
                # Convert roughness from millimeters to meters.
                pipe["roughness"] = 0.001 * parse(Float64, current[6])
            elseif data["head_loss"] == "H-W" # If Hazen-Williams head loss is used...
                # Retain the original value (unitless).
                pipe["roughness"] = parse(Float64, current[6])
            end
        elseif data["flow_units"] == "GPM" # If gallons per minute...
            # Convert length from feet to meters.
            pipe["length"] = 0.3048 * parse(Float64, current[4])

            # Convert diameter from inches to meters.
            pipe["diameter"] = 0.0254 * parse(Float64, current[5])

            if data["head_loss"] == "D-W" # If Darcy-Weisbach head loss is used...
                # Convert roughness from millifeet to meters.
                pipe["roughness"] = 3.048e-4 * parse(Float64, current[6])
            elseif data["head_loss"] == "H-W" # If Hazen-Williams head loss is used...
                # Retain the original value (unitless).
                pipe["roughness"] = parse(Float64, current[6])
            end
        else
            error("Could not find a valid \"units\" option type.")
        end

        # Parse minor loss data (unitless).
        pipe["minor_loss"] = parse(Float64, current[7])

        # Derive important metadata from existing data.
        pipe["has_valve"] = uppercase(current[8]) in ["CV", "CLOSED"]
        pipe["flow_direction"] = uppercase(current[8]) == "CV" ? FLOW_DIRECTION_POSITIVE : FLOW_DIRECTION_UNKNOWN

        # Add a temporary index to be used in the data dictionary.
        pipe["index"] = string(index += 1)

        # Append the pipe entry to the data dictionary.
        data["pipe"][current[1]] = pipe
    end

    # Replace with new dictionaries that use integer component indices.
    data["pipe"] = _transform_component_indices(data["pipe"])
end

function _read_pump!(data::Dict{String,<:Any})
    # Initialize dictionaries associated with pumps.
    data["pump"] = Dict{String,Dict{String,Any}}()
    data["time_series"]["pump"] = Dict{String,Any}()

    # Initialize a temporary index to be updated while parsing.
    index::Int = 0

    # Loop over all lines in the [PUMPS] section and parse each.
    for (line_number, line) in data["section"]["[PUMPS]"]
        current = split(split(line, ";")[1])
        length(current) == 0 && continue # Skip.

        # Initialize the pipe entry to be added.
        pump = Dict{String,Any}("name" => current[1], "status" => STATUS_UNKNOWN)
        pump["source_id"] = ["pump", current[1]]
        pump["node_fr"] = data["node_map"][current[2]]
        pump["node_to"] = data["node_map"][current[3]]

        # Loop over remaining entries and store remaining properties.
        for i in range(4; stop = length(current), step = 2)
            if uppercase(current[i]) == "POWER"
                # Memento.error(_LOGGER, "Constant power pumps are not supported.")
                pump["pump_type"] = PUMP_CONSTANT_POWER

                # Parse the head of the reservoir (in meters).
                if data["flow_units"] == "LPS" || data["flow_units"] == "CMH"
                    # Conver power from kilowatts to watts.
                    pump["power"] = 1.0e3 * parse(Float64, current[i+1])
                elseif data["flow_units"] == "GPM" # If gallons per minute...
                    # Convert power from horsepower to watts.
                    pump["power"] = 745.7 * parse(Float64, current[i+1])
                else
                    Memento.error(_LOGGER, "Could not find a valid \"units\" option type.")
                end
            elseif uppercase(current[i]) == "HEAD"
                # Obtain and scale head-versus-flow pump curve.
                flow = first.(data["curve"][current[i+1]]) # Flow.
                head = last.(data["curve"][current[i+1]]) # Head.

                if data["flow_units"] == "LPS" # If liters per second...
                    # Convert from liters per second to cubic meters per second.
                    flow *= 1.0e-3
                elseif data["flow_units"] == "CMH" # If cubic meters per hour...
                    # Convert from cubic meters per hour to cubic meters per second.
                    flow *= inv(3600.0)
                elseif data["flow_units"] == "GPM" # If gallons per minute...
                    # Convert from gallons per minute to cubic meters per second.
                    flow *= 6.30902e-5

                    # Convert head from feet to meters.
                    head *= 0.3048
                else
                    Memento.error(_LOGGER, "Could not find a valid \"units\" option type.")
                end

                # Curve of head (meters) versus flow (cubic meters per second).
                pump["head_curve"] = Array([(flow[j], head[j]) for j = 1:length(flow)])     
            elseif uppercase(current[i]) != "HEAD"
                Memento.error(_LOGGER, "Pump keyword in INP file not recognized.")
            end
        end

        # Flow is always in the positive direction for pumps.
        pump["flow_direction"] = FLOW_DIRECTION_POSITIVE

        # Add a temporary index to be used in the data dictionary.
        pump["index"] = string(index += 1)

        # Append the pump entry to the data dictionary.
        data["pump"][current[1]] = pump
    end

    # Replace with new dictionaries that use integer component indices.
    data["pump"] = _transform_component_indices(data["pump"])
    data["time_series"]["pump"] = _transform_time_series_indices(
        data["pump"], data["time_series"]["pump"])
end

function _read_reservoir!(data::Dict{String,<:Any})
    # Initialize dictionaries associated with reservoirs.
    data["reservoir"] = Dict{String,Dict{String,Any}}()
    data["time_series"]["reservoir"] = Dict{String,Any}()

    # Initialize a temporary index to be updated while parsing.
    index::Int = 0

    # Loop over all lines in the [RESERVOIRS] section and parse each.
    for (line_number, line) in data["section"]["[RESERVOIRS]"]
        current = split(split(line, ";")[1])
        length(current) == 0 && continue # Skip.

        # Initialize the reservoir entry to be added.
        reservoir = Dict{String,Any}("name" => current[1], "status" => STATUS_ACTIVE)
        reservoir["source_id"] = ["reservoir", current[1]]
        reservoir["dispatchable"] = false

        # Parse the head of the reservoir (in meters).
        if data["flow_units"] == "LPS" || data["flow_units"] == "CMH"
            # Retain the original value (in meters).
            reservoir["head_nominal"] = parse(Float64, current[2])
        elseif data["flow_units"] == "GPM" # If gallons per minute...
            # Convert head from feet to meters.
            reservoir["head_nominal"] = 0.3048 * parse(Float64, current[2])
        else
            Memento.error(_LOGGER, "Could not find a valid \"units\" option type.")
        end

        # Parse the name of the pattern used to scale head at the reservoir.
        pattern = length(current) > 2 ? current[3] : nothing

        # Scale the parsed head by the specified pattern, if it exists. Also, assume the
        # reservoir's elevation is equal to the minimum value of head presented in the data.
        if pattern !== nothing && length(data["pattern"][pattern]) > 1
            head = reservoir["head_nominal"] .* data["pattern"][pattern]
            push!(head, maximum(head))

            entry = Dict{String,Array{Float64}}("head_nominal" => head)
            entry["head_min"] = entry["head_nominal"]
            entry["head_max"] = entry["head_nominal"]
            data["time_series"]["reservoir"][current[1]] = entry

            reservoir["head_min"] = minimum(entry["head_min"])
            reservoir["head_max"] = maximum(entry["head_max"])
            reservoir["elevation"] = minimum(head)
        elseif pattern !== nothing && pattern != "1"
            reservoir["head_nominal"] *= data["pattern"][pattern][1]
            reservoir["head_min"] = reservoir["head_nominal"]
            reservoir["head_max"] = reservoir["head_nominal"]
            reservoir["elevation"] = reservoir["head_nominal"]
        else
            reservoir["elevation"] = reservoir["head_nominal"]
            reservoir["head_min"] = reservoir["head_nominal"]
            reservoir["head_max"] = reservoir["head_nominal"]
        end

        # Add a temporary index to be used in the data dictionary.
        reservoir["index"] = string(index += 1)

        # Set the unit cost for extraction and treatment for outgoing flow from a reservoir
        # (in units of currency per cubic meter per second). (The default cost is zero.)
        reservoir["flow_cost"] = 0.0 # EPANET does not provide this data.

        # Append the reservoir entry to the data dictionary.
        data["reservoir"][current[1]] = reservoir
    end

    # Replace with new dictionaries that use integer component indices.
    data["reservoir"] = _transform_component_indices(data["reservoir"])
    #data["time_series"]["reservoir"] = _transform_time_series_indices(
    #    data["reservoir"], data["time_series"]["reservoir"])
end

function _read_tank!(data::Dict{String,<:Any})
    # Initialize dictionaries associated with tanks.
    data["tank"] = Dict{String,Dict{String,Any}}()

    # Initialize a temporary index to be updated while parsing.
    index::Int = 0

    # Loop over all lines in the [TANKS] section and parse each.
    for (line_number, line) in data["section"]["[TANKS]"]
        current = split(split(line, ";")[1])
        length(current) == 0 && continue # Skip.

        # Initialize the tank entry to be added.
        tank = Dict{String,Any}("name" => current[1], "status" => STATUS_ACTIVE)
        tank["source_id"] = ["tank", current[1]]
        tank["dispatchable"] = false

        # Store all measurements associated with tanks in metric units.
        if data["flow_units"] == "LPS" || data["flow_units"] == "CMH"
            # Retain the original length values (in meters).
            tank["elevation"] = parse(Float64, current[2])
            tank["init_level"] = parse(Float64, current[3])
            tank["min_level"] = parse(Float64, current[4])
            tank["max_level"] = parse(Float64, current[5])
            tank["diameter"] = parse(Float64, current[6])

            # Retain the original minimum volume (in cubic meters).
            tank["min_vol"] = parse(Float64, current[7])
        elseif data["flow_units"] == "GPM" # If gallons per minute...
            # Convert length values from feet to meters.
            tank["elevation"] = 0.3048 * parse(Float64, current[2])
            tank["init_level"] = 0.3048 * parse(Float64, current[3])
            tank["min_level"] = 0.3048 * parse(Float64, current[4])
            tank["max_level"] = 0.3048 * parse(Float64, current[5])
            tank["diameter"] = 0.3048 * parse(Float64, current[6])

            # Convert minimum volume from cubic feet to cubic meters.
            tank["min_vol"] = 0.3048^3 * parse(Float64, current[7])
        else
            Memento.error(_LOGGER, "Could not find a valid \"Units\" option type.")
        end

        # Add a temporary index to be used in the data dictionary.
        tank["index"] = string(index += 1)

        # Append the tank entry to the data dictionary.
        data["tank"][current[1]] = tank
    end

    # Replace with new dictionaries that use integer component indices.
    data["tank"] = _transform_component_indices(data["tank"])
end

function _read_time!(data::Dict{String,<:Any})
    # Loop over all lines in the [TIMES] section and parse each.
    for (line_number, line) in data["section"]["[TIMES]"]
        current = split(split(line, ";")[1])
        length(current) == 0 && continue # Skip.

        if uppercase(current[1]) == "DURATION"
            data["duration"] = _string_time_to_seconds(current[2])
        elseif uppercase(current[1]) == "HYDRAULIC" && uppercase(current[2]) == "TIMESTEP"
            data["hydraulic_time_step"] = _string_time_to_seconds(current[3])
        elseif uppercase(current[1]) == "PATTERN" && uppercase(current[2]) == "TIMESTEP"
            data["pattern_time_step"] = _string_time_to_seconds(current[3])
        elseif uppercase(current[1]) == "CLOCKTIME"
            time_format = length(current) > 3 ? uppercase(current[4]) : "AM"
            data["start_time"] = _clock_time_to_seconds(current[3], time_format)
        end
    end
end

function _read_title!(data::Dict{String,<:Any})
    lines = Array{String}[]

    for (line_number, line) in data["section"]["[TITLE]"]
        line = split(line, ";")[1]
        current = split(line)

        if length(current) == 0
            continue
        else
            lines = vcat(lines, line)
        end
    end

    return data["name"] = join(lines, " ")
end

function _read_regulator!(data::Dict{String,<:Any})
    # Initialize dictionaries associated with valves.
    data["regulator"] = Dict{String,Dict{String,Any}}()

    # Initialize a temporary index to be updated while parsing.
    index::Int = 0

    # Loop over all lines in the [VALVES] section and parse each.
    for (line_number, line) in data["section"]["[VALVES]"]
        current = split(split(line, ";")[1])
        length(current) == 0 && continue # Skip.

        # Initialize the valve entry to be added.
        valve = Dict{String,Any}("name" => current[1], "status" => STATUS_UNKNOWN)
        valve["node_fr"] = data["node_map"][current[2]]
        valve["node_to"] = data["node_map"][current[3]]

        # Parse the valve type and throw an error if not supported.
        if uppercase(current[5]) == "PRV"
            valve["source_id"] = ["regulator", current[1]]
            valve["flow_direction"] = FLOW_DIRECTION_POSITIVE
        elseif uppercase(current[5]) == "TCV"
            valve["source_id"] = ["throttle_control_valve", current[1]]
            Memento.error(_LOGGER, "Valves of type $(current[5]) are not supported.")
        else
            Memento.error(_LOGGER, "Valves of type $(current[5]) are not supported.")
        end

        # Store all measurements associated with pipes in metric units.
        if data["flow_units"] == "LPS" || data["flow_units"] == "CMH"
            # Retain the original setting value (in meters).
            valve["setting"] = parse(Float64, current[6])

            # Convert diameter from millimeters to meters.
            valve["diameter"] = 0.001 * parse(Float64, current[4])
        elseif data["flow_units"] == "GPM" # If gallons per minute...
            # Convert setting from PSI to meters.
            valve["setting"] = inv(1.421970206324753) * parse(Float64, current[6])

            # Convert diameter from inches to meters.
            valve["diameter"] = 0.0254 * parse(Float64, current[4])
        else
            Memento.error(_LOGGER, "Could not find a valid \"units\" option type.")
        end

        # Parse minor loss data (unitless).
        valve["minor_loss"] = parse(Float64, current[7])

        # Add a temporary index to be used in the data dictionary.
        valve["index"] = string(index += 1)

        # Append the valve entry to the data dictionary.
        data[valve["source_id"][1]][current[1]] = valve
    end

    # Replace with a new dictionary that uses integer component indices.
    data["regulator"] = _transform_component_indices(data["regulator"])
end

"""
    epanet_to_watermodels!(epanet_data; import_all=false)

Converts data parsed from an EPANET file, passed by `epanet_data` into a format suitable for
internal WaterModels use. Imports all data from the EPANET file if `import_all` is true.
"""
function epanet_to_watermodels!(data::Dict{String,<:Any}; import_all::Bool = false)
    drop_zero_demands!(data) # Drop demands of zero from nodes.
    # convert_short_pipes!(data) # Convert pipes that are short to short pipes and valves.
    # _add_valves_to_tanks!(data) # Ensure that shutoff valves are connected to tanks.
    add_valves_from_pipes!(data) # Convert pipes with valves to pipes *and* valves.
    drop_pipe_flags!(data) # Drop irrelevant pipe attributes.
end


function _get_max_comp_id(data::Dict{String,<:Any}, comp::String)
    return length(data[comp]) > 0 ? maximum([x["index"] for (i, x) in data[comp]]) : 0
end


function _get_max_link_id(data::Dict{String,<:Any})
    return maximum(maximum([x["index"] for (i, x) in data[t]]) for t in _LINK_COMPONENTS)
end


function _get_max_node_id(data::Dict{String,<:Any})
    return maximum([x["index"] for (i, x) in data["node"]])
end


function convert_short_pipes!(data::Dict{String, <:Any})
    wm_data = get_wm_data(data)
    head_loss, viscosity = wm_data["head_loss"], wm_data["viscosity"]
    func! = x -> _convert_short_pipes!(x, head_loss, viscosity)
    apply_wm!(func!, data; apply_to_subnetworks = true)
end


function _convert_short_pipes!(data::Dict{String,<:Any}, head_loss::String, viscosity::Float64)
    exponent = _get_exponent_from_head_loss_form(head_loss)
    max_flow_exp = abs(calc_capacity_max(data))^exponent

    wm_data = get_wm_data(data)
    head_transform = _calc_head_per_unit_transform(wm_data)
    base_length = 1.0 / _calc_length_per_unit_transform(wm_data)(1.0)
    base_mass = 1.0 / _calc_mass_per_unit_transform(wm_data)(1.0)
    base_time = 1.0 / _calc_time_per_unit_transform(wm_data)(1.0)

    for (a, pipe) in data["pipe"]
        r = _calc_pipe_resistance(pipe, head_loss, viscosity, base_length, base_mass, base_time)
        dh_max = pipe["length"] * r * max_flow_exp

        # If maximum head loss is less than one centimeter...
        if dh_max <= head_transform(0.01)
            # Delete unnecessary fields.
            delete!(pipe, "diameter")
            delete!(pipe, "length")
            delete!(pipe, "roughness")

            if haskey(pipe, "has_valve") && pipe["has_valve"]
                # Transform the pipe into a valve.
                delete!(pipe, "has_valve")
                data["valve"][a] = deepcopy(pipe)
            else
                # Transform the pipe into a short pipe.
                delete!(pipe, "has_valve")
                data["short_pipe"][a] = deepcopy(pipe)
            end

            # Delete the original pipe component.
            delete!(data["pipe"], a)
        end
    end
end


function _has_zero_demand(demand::Dict{String, <:Any})
    zero_flow_min = demand["flow_min"] == 0.0
    zero_flow_max = demand["flow_max"] == 0.0
    return zero_flow_min && zero_flow_max
end


function drop_zero_demands!(data::Dict{String, <:Any})
    apply_wm!(_drop_zero_demands!, data; apply_to_subnetworks = true)
end


function _drop_zero_demands!(data::Dict{String,<:Any})
    for (i, demand) in filter(x -> _has_zero_demand(x.second), data["demand"])
        delete!(data["demand"], i)

        if haskey(data, "time_series") && haskey(data["time_series"], "demand")
            delete!(data["time_series"]["demand"], i)
        end
    end
end


function drop_pipe_flags!(data::Dict{String, <:Any})
    apply_wm!(_drop_pipe_flags!, data; apply_to_subnetworks = true)
end


function _drop_pipe_flags!(data::Dict{String,<:Any})
    # Drop valve statuses to parsed pipes.
    for (a, pipe) in data["pipe"]
        delete!(pipe, "has_valve")
    end
end


function _add_valves_to_tanks!(data::Dict{String,<:Any})
    # Modify the network for standard modeling of tanks.
    for (i, tank) in data["tank"]
        # Get the original index of the tank node.
        original_tank_node = tank["node"]

        # Create a new dummy node at which the tank will now be placed.
        dummy_node_id = _get_max_comp_id(data, "node") + 1
        dummy_node = deepcopy(data["node"][string(tank["node"])])
        dummy_node["index"], dummy_node["name"] = dummy_node_id, string(dummy_node_id)
        data["node"][string(dummy_node_id)] = dummy_node
        tank["node"] = dummy_node_id # The tank now resides at the dummy node.

        # Add a valve that connects the dummy (tank) node to the original node.
        v_id = _get_max_comp_id(data, "valve") + 1
        valve = Dict{String,Any}("name" => string(v_id), "status" => STATUS_UNKNOWN, "index" => v_id)
        valve["source_id"] = AbstractString["valve", string(v_id)]
        valve["node_fr"], valve["node_to"] = original_tank_node, dummy_node_id
        valve["flow_direction"], valve["minor_loss"] = FLOW_DIRECTION_UNKNOWN, 0.0

        # Add the valve to the data object.
        data["valve"][string(v_id)] = valve
    end
end


function add_valves_from_pipes!(data::Dict{String, <:Any})
    apply_wm!(_add_valves_from_pipes!, data; apply_to_subnetworks = true)
end

function _add_valves_from_pipes!(data::Dict{String,<:Any})
    # Add valves to parsed pipes that have the attribute "has_valve" equal to true.
    for (a, pipe) in filter(x -> get(x.second, "has_valve", false), data["pipe"])
        # Create a dummy node to which the valve and pipe will be attached.
        dummy_node_id = _get_max_comp_id(data, "node") + 1
        dummy_node = deepcopy(data["node"][string(pipe["node_fr"])])
        dummy_node["index"], dummy_node["name"] = dummy_node_id, string(dummy_node_id)

        # Create the valve that will connect to the dummy node.
        valve = Dict{String,Any}("index" => pipe["index"], "name" => pipe["name"])
        valve["source_id"], valve["status"] = pipe["source_id"], STATUS_UNKNOWN
        valve["flow_direction"] = pipe["flow_direction"]
        valve["minor_loss"] = pipe["minor_loss"]

        # Modify the start and ends nodes of the valve and pipe.
        valve["node_fr"], valve["node_to"] = pipe["node_fr"], dummy_node_id
        pipe["node_fr"], pipe["node_to"] = dummy_node_id, pipe["node_to"]

        # Store the valve object in the valve dictionary.
        data["valve"][a], data["node"][string(dummy_node_id)] = valve, dummy_node
    end
end


_INP_SECTIONS = [
    "[OPTIONS]",
    "[TITLE]",
    "[JUNCTIONS]",
    "[RESERVOIRS]",
    "[TANKS]",
    "[PIPES]",
    "[PUMPS]",
    "[VALVES]",
    "[EMITTERS]",
    "[CURVES]",
    "[PATTERNS]",
    "[ENERGY]",
    "[STATUS]",
    "[CONTROLS]",
    "[RULES]",
    "[DEMANDS]",
    "[QUALITY]",
    "[REACTIONS]",
    "[SOURCES]",
    "[MIXING]",
    "[TIMES]",
    "[REPORT]",
    "[COORDINATES]",
    "[VERTICES]",
    "[LABELS]",
    "[BACKDROP]",
    "[TAGS]",
]
