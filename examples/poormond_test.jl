import Xpress
import JuMP
using Revise
using WaterModels


modification_path = "/Users/dsigler/PowerWaterExamples/data/water/benchmarks/Poormond/modifications.json"
network_path = "/Users/dsigler/PowerWaterExamples/data/water/benchmarks/Poormond/Poormond-12_Steps-Day_21.inp"
network = WaterModels.parse_file(network_path; skip_correct = true)
modifications = WaterModels.parse_file(modification_path; skip_correct = true)
WaterModels._IM.update_data!(network, modifications)
WaterModels.correct_network_data!(network)
network_mn = WaterModels.make_multinetwork(network);
ext = Dict{Symbol, Any}(:pipe_breakpoints => 3, :pump_breakpoints => 3)
wm = WaterModels.instantiate_model(network_mn , PWLRDWaterModel, build_mn_owf; ext = ext)
#relax_all_binary_variables!(wm)
xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>60.0);
result_LP = WaterModels.optimize_model!(wm, optimizer = xp)





