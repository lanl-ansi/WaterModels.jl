#import Xpress
import Gurobi
import JuMP
using Revise
using WaterModels
 
network = parse_file("/scratch/dsigler/PowerWaterExamples/data/water/6-10_pumps/Richmond_standard.inp");
network_mn = WaterModels.make_multinetwork(network);
ext = Dict{Symbol, Any}(:pipe_breakpoints => 5, :pump_breakpoints => 5)
wm = WaterModels.instantiate_model(network_mn , PWLRDWaterModel, build_mn_owf; ext = ext)
#relax_all_binary_variables!(wm)

env = Gurobi.Env();
gurobi = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "TimeLimit" => 60.0);
result= WaterModels.optimize_model!(wm, optimizer = gurobi);