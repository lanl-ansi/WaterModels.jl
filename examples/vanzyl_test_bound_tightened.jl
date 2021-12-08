using Pkg
Pkg.activate("/Users/dsigler/WaterModels.jl")

import Xpress
import JuMP
using Revise
using WaterModels

# network_path = "/Users/btasseff/Desktop/Van_Zyl-BT.json";
# network_mn = WaterModels.parse_file(network_path; skip_correct = true);
# WaterModels.set_flow_partitions_si!(network_mn, 1.0, 1.0e-4);
# result = WaterModels.solve_mn_owf(network_mn, WaterModels.PWLRDWaterModel, gurobi);

#relax_all_binary_variables!(wm)
xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>300.0);
network_path = "./examples/data/json/Van_Zyl-BT.json";
network_mn = WaterModels.parse_file(network_path; skip_correct = true);
WaterModels.set_flow_partitions_si!(network_mn, 5.0, 1.0e-4);
wm = WaterModels.instantiate_model(network_mn, WaterModels.PWLRDWaterModel, WaterModels.build_mn_owf);
result = WaterModels.optimize_model!(wm, optimizer = xp);

pos_flow_array = []
neg_flow_array =[]
# h_4_array = []
# h_13_array = []
for time in  1:4
#     println("")
#     println("time ",time)
#     h_4 = JuMP.value(wm.var[:it][:wm][:nw][time][:h][4])
#     h_13 = JuMP.value(wm.var[:it][:wm][:nw][time][:h][13])
    pos_flow_4_13 = JuMP.value(wm.var[:it][:wm][:nw][time][:qp_pipe][2])
    neg_flow_4_13 = JuMP.value(wm.var[:it][:wm][:nw][time][:qn_pipe][2])
    push!(pos_flow_array,pos_flow_4_13)
    push!(neg_flow_array,neg_flow_4_13)
#     push!(h_4_array,h_4)
#     push!(h_13_array,h_13)
#     @show h_4
#     @show h_13
#     @show pos_flow_4_13
#     @show neg_flow_4_13
#     println("")
end
@show pos_flow_array
@show neg_flow_array
# @show h_4_array
# @show h_13_array


