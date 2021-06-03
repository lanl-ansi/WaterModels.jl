import Gurobi
import JuMP

using Revise
using WaterModels
 
network = parse_file("examples/data/epanet/van_zyl.inp");
network_mn = WaterModels.make_multinetwork(network);

ext = Dict{Symbol, Any}(:pipe_breakpoints => 5, :pump_breakpoints => 5);
wm = WaterModels.instantiate_model(network_mn , PWLRDWaterModel, build_mn_owf; ext = ext);
relax_all_binary_variables!(wm);

env = Gurobi.Env();
gurobi = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "TimeLimit" => 15.0);
result_lp = WaterModels.optimize_model!(wm, optimizer = gurobi);

split_networks = split_multinetwork(network_mn,
    [["1", "2", "3", "4", "5", "6"],
     ["6", "7", "8","9", "10", "11", "12"],
     ["12", "13", "14", "15", "16", "17"],
     ["17", "18", "19", "20", "21", "22", "23", "24"]]);


# Initialize empty array of WaterModels objects.
wms = Array{AbstractWaterModel}([]) 
count = 1

for subnetwork in split_networks
    global count

    if count == 1
        wm = WaterModels.instantiate_model(subnetwork, PWLRDWaterModel, build_mn_owf_part_start_int; ext = ext)
        append!(wms, [wm])
    else
        wm = WaterModels.instantiate_model(subnetwork, PWLRDWaterModel, build_mn_owf_part_int; ext = ext)
        append!(wms, [wm])
    end

    count += 1
end

for part in wms
    println("\n\n\n\n\n\nstarting new problem\n\n\n\n\n\n")
    result = WaterModels.optimize_model!(part; optimizer = gurobi)
end

