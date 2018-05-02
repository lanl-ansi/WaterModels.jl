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

    return epanet_dict
end
