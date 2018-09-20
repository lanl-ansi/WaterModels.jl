import JSON

function parse_file(file::String)
    if endswith(file, ".inp")
        network_data = WaterModels.parse_epanet_file(file)
    elseif endswith(file, ".json")
        network_data = JSON.parsefile(file)
        network_data["per_unit"] = false
    else
        error("'" + file + "' is not a valid file type.")
    end
    
    return network_data
end

function print_solution(solution::Dict)
    InfrastructureModels.print_summary(solution)
end
