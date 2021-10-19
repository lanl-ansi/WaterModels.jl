import Xpress
import JuMP
using Revise
using WaterModels
 
network = parse_file("examples/data/epanet/van_zyl.inp");
network_mn = WaterModels.make_multinetwork(network);
ext = Dict{Symbol, Any}(:pipe_breakpoints => 5, :pump_breakpoints => 5)
wm = WaterModels.instantiate_model(network_mn , PWLRDWaterModel, build_mn_owf; ext = ext)
#relax_all_binary_variables!(wm)
xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>3600.0);
result_LP = WaterModels.optimize_model!(wm, optimizer = xp)

# split_networks = split_multinetwork(network_mn, 
# [["1", "2", "3", "4", "5","6","7"], 
# ["7", "8","9","10","11","12","13"],
# ["13","14","15","16","17","18"],
# ["18","19","20","21","22","23","24"]]);

 
# wms = Array{AbstractWaterModel}([]) # Initialize empty array of WaterModels objects.
# count = 1

# for subnetwork in split_networks
#     global count
#     if count == 1
#         ext = Dict{Symbol, Any}(:pipe_breakpoints => 5, :pump_breakpoints => 5)
#         wm = WaterModels.instantiate_model(subnetwork, PWLRDWaterModel, build_mn_owf_part_start_v2; ext = ext)
#         append!(wms, [wm])
#     elseif count == 4 
#         ext = Dict{Symbol, Any}(:pipe_breakpoints => 5, :pump_breakpoints => 5)
#         wm = WaterModels.instantiate_model(subnetwork, PWLRDWaterModel, build_mn_owf_part_end_v2; ext = ext)
#         append!(wms, [wm])

#     else
#         ext = Dict{Symbol, Any}(:pipe_breakpoints => 5, :pump_breakpoints => 5)
#         wm = WaterModels.instantiate_model(subnetwork, PWLRDWaterModel, build_mn_owf_part_middle_v2; ext = ext)
#         append!(wms, [wm])
#     end
#     count += 1
# end

# xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>300.0);

# for part in wms

#     println("\n\n\n\n\n\nstarting new problem\n\n\n\n\n\n")

#     xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>15.0);
 
#     # JuMP.set_optimizer(part.model, Xpress.Optimizer)
#     # JuMP.set_optimizer_attributes(part.model, "MAXTIME"=>300.0)
#     # JuMP.optimize!(part.model)

#     result = WaterModels.optimize_model!(part, optimizer = xp)
# end

