#Poormond-12_Steps-Day_21-BT-TS.json

using Pkg
Pkg.activate("/Users/dsigler/WaterModels.jl")

import Xpress
import JuMP
using Revise
using WaterModels




#relax_all_binary_variables!(wm)
xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>30.0);
network_path = "./examples/data/json/Poormond-12_Steps-Day_21-BT-TS.json";
network_mn = WaterModels.parse_file(network_path; skip_correct = true);
WaterModels.set_flow_partitions_si!(network_mn, 5.0, 1.0e-4);
wm = WaterModels.instantiate_model(network_mn, WaterModels.LRDWaterModel, WaterModels.build_mn_owf);
WaterModels.optimize_model!(wm, optimizer = xp);




