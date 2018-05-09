############################################################################
#                                                                          #
# This file provides functionality for interfacing with EPANET data files. #
# See https://github.com/OpenWaterAnalytics/EPANET/wiki/Input-File-Format  #
# for a thorough description of the format and its components.             #
#                                                                          #
############################################################################

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

    # Correct the type of the title entry within the dictionary.
    if length(keys(epanet_dict["title"])) > 0
        epanet_dict["title"] = first(keys(epanet_dict["title"]))
    else
        epanet_dict["title"] = ""
    end

    # Parse relevant data into a more structured format.
    epanet_dict["junctions"] = parse_junctions(epanet_dict["junctions"])

    epanet_dict["multinetwork"] = false

    return epanet_dict
end

function allequal(x) 
    return all(y->y == x[1], x)
end

function parse_junctions(data::Dict{String, Array})
    # Specify the data types for the junction data.
    columns = Dict("demand" => Float64, "elev" => Float64,
                   "id" => String, "pattern" => String)

    # Ensure the arrays describing junction data are all of equal lengths.
    @assert(allequal([length(data[column]) for column in keys(columns)]))

    # Build an array of junction dictionaries with the correct data types.
    num_junctions = length(data["id"])
    #junctions = [Dict(c => parse(v, data[c][i]) for (c, v) in columns) for i = 1:num_junctions]
    #println(junctions)

    return data
end
