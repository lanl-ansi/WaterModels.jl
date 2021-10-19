using Pkg
Pkg.activate("/Users/dsigler/WaterModels.jl")

import Xpress
import JuMP
using Revise
using WaterModels
 
#network = parse_file("/scratch/dsigler/PowerWaterExamples/data/water/6-10_pumps/Richmond_skeleton_raised.inp");
network = parse_file("./examples/data/epanet/van_zyl.inp");
network_mn = WaterModels.make_multinetwork(network);
set_flow_partitions_si!(network_mn, 10.0, 1.0e-4)
#ext = Dict{Symbol, Any}(:pipe_breakpoints => 5, :pump_breakpoints => 5)
#wm = WaterModels.instantiate_model(network_mn , PWLRDWaterModel, build_mn_owf; ext = ext)
wm = WaterModels.instantiate_model(network_mn , PWLRDWaterModel, build_mn_owf;)
#relax_all_binary_variables!(wm)
xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>30.0);
result_LP = WaterModels.optimize_model!(wm, optimizer = xp)