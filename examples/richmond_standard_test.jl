import Xpress
import JuMP
using Revise
using WaterModels
 
network = parse_file("/Users/dsigler/PowerWaterExamples/data/water/6-10_pumps/Richmond_standard.inp");
network_mn = WaterModels.make_multinetwork(network);
ext = Dict{Symbol, Any}(:pipe_breakpoints => 5, :pump_breakpoints => 5)
wm = WaterModels.instantiate_model(network_mn , PWLRDWaterModel, build_mn_owf; ext = ext)
#relax_all_binary_variables!(wm)
xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>3600.0);
result_LP = WaterModels.optimize_model!(wm, optimizer = xp)

