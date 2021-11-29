import Gurobi
#import Xpress
import JuMP
using Revise
using WaterModels


modification_path = "/scratch/dsigler/PowerWaterExamples/data/water/benchmarks/Poormond/modifications.json"
network_path = "/scratch/dsigler/PowerWaterExamples/data/water/benchmarks/Poormond/Poormond-24_Steps-Day_21.inp"
network = WaterModels.parse_file(network_path; skip_correct = true)
modifications = WaterModels.parse_file(modification_path; skip_correct = true)
WaterModels._IM.update_data!(network, modifications)
WaterModels.correct_network_data!(network)
network_mn = WaterModels.make_multinetwork(network);
ext = Dict{Symbol, Any}(:pipe_breakpoints => 3, :pump_breakpoints => 3)
wm = WaterModels.instantiate_model(network_mn , PWLRDWaterModel, build_mn_owf; ext = ext)
env = Gurobi.Env();
gurobi = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "TimeLimit" => 3600.0);
result= WaterModels.optimize_model!(wm, optimizer = gurobi);





