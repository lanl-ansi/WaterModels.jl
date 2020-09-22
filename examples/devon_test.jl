
import Xpress
import JuMP
using Revise
using WaterModels

xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MIPABSSTOP"=>0.0);
network = parse_file("examples/data/epanet/van_zyl.inp");
ext = Dict(:pipe_breakpoints=>10, :pump_breakpoints=>5);
#run_obbt_owf!(network, xp, model_type = LRDWaterModel, time_limit=5.0, ext=ext, relaxed=false);
network_mn = WaterModels.make_multinetwork(network);
make_tank_start_dispatchable!(network_mn);
xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>60.0);
result = run_mn_owf(network_mn, LRDWaterModel, xp; ext=ext);
