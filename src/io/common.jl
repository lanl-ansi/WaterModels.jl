function parse_file(file::String)
    if endswith(file, ".inp")
        network_data = WaterModels.parse_epanet_file(file)
    end
    
    return network_data
end

function print_solution(solution::Dict)
    InfrastructureModels.print_summary(solution)
end
