"""
Implements an algorithm to determine the globally optimal solution to a network
design problem. (This is a reproduction of Algorithm 1 in Raghunathan (2013).)
"""
function solve_global(network_path::String, problem_path::String)
    network = WaterModels.parse_file(network_path)
    modifications = WaterModels.parse_file(problem_path)
    InfrastructureModels.update_data!(network, modifications)
end
