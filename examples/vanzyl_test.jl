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
# set_points = [0.1519191109426697,0.12815817888078826,0.1661881436210814,0.18622668547036325,0.10986057366303964,0.1750347111100875]
# for time in [6]#1:(length(set_points)-4) # 1 ok 3 ok 4 ok
#     vref = wm.var[:it][:wm][:nw][time][:qp_pipe][2]
#     JuMP.fix(vref, set_points[time]; force = true)
# end
result_LP = WaterModels.optimize_model!(wm, optimizer = xp)
pos_flow_array = []
neg_flow_array =[]
h_4_array = []
h_13_array = []
for time in  1:4
    println("")
    println("time ",time)
    h_4 = JuMP.value(wm.var[:it][:wm][:nw][time][:h][4])
    h_13 = JuMP.value(wm.var[:it][:wm][:nw][time][:h][13])
    pos_flow_4_13 = JuMP.value(wm.var[:it][:wm][:nw][time][:qp_pipe][2])
    neg_flow_4_13 = JuMP.value(wm.var[:it][:wm][:nw][time][:qn_pipe][2])
    push!(pos_flow_array,pos_flow_4_13)
    push!(neg_flow_array,neg_flow_4_13)
    push!(h_4_array,h_4)
    push!(h_13_array,h_13)
    @show h_4
    @show h_13
    @show pos_flow_4_13
    @show neg_flow_4_13
    println("")
end
@show pos_flow_array
@show neg_flow_array
@show h_4_array
@show h_13_array


time 1
h_4 = 1.0247915558722451
h_13 = 0.8793539521371608
pos_flow_4_13 = 0.2411775428041607
neg_flow_4_13 = 0.0


time 2
h_4 = 0.9894384162847942
h_13 = 0.8225119671669396
pos_flow_4_13 = 0.2598190603862963
neg_flow_4_13 = 0.0


time 3
h_4 = 0.8635675455329749
h_13 = 0.7697081113168751
pos_flow_4_13 = 0.19643361464907033
neg_flow_4_13 = 0.0


time 4
h_4 = 0.7681818181818182
h_13 = 0.7681818181818182
pos_flow_4_13 = 0.0
neg_flow_4_13 = 0.0

pos_flow_array = Any[0.2411775428041607, 0.2598190603862963, 0.19643361464907033, 0.0]
neg_flow_array = Any[0.0, 0.0, 0.0, 0.0]
h_4_array = Any[1.0247915558722451, 0.9894384162847942, 0.8635675455329749, 0.7681818181818182]
h_13_array = Any[0.8793539521371608, 0.8225119671669396, 0.7697081113168751, 0.7681818181818182]
4-element Vector{Any}:
 0.8793539521371608
 0.8225119671669396
 0.7697081113168751
 0.7681818181818182
