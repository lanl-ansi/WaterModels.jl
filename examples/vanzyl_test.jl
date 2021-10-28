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

for time in  1:6
    println("")
    println("time ",time)
    h_4 = JuMP.value(wm.var[:it][:wm][:nw][time][:h][4])
    h_13 = JuMP.value(wm.var[:it][:wm][:nw][time][:h][13])
    pos_flow_4_13 = JuMP.value(wm.var[:it][:wm][:nw][time][:qp_pipe][2])
    neg_flow_4_13 = JuMP.value(wm.var[:it][:wm][:nw][time][:qn_pipe][2])

    @show h_4
    @show h_13
    @show pos_flow_4_13
    @show neg_flow_4_13
    println("")
end


# time 1
# h_4 = 0.7681818181818182
# h_13 = 0.7681818181818182
# pos_flow_4_13 = 0.0
# neg_flow_4_13 = 0.0


# time 2
# h_4 = 1.0224107797360427
# h_13 = 0.8718787201021464
# pos_flow_4_13 = 0.24559697007060888
# neg_flow_4_13 = 0.0


# time 3
# h_4 = 1.0014814907696807
# h_13 = 0.8156485165421574
# pos_flow_4_13 = 0.2686054280077992
# neg_flow_4_13 = 0.0


# time 4
# h_4 = 0.8604757199477432
# h_13 = 0.7651654378411213
# pos_flow_4_13 = 0.19769222142753703
# neg_flow_4_13 = 0.0


# time 5
# h_4 = 0.8607492101967075
# h_13 = 0.7655672641620435
# pos_flow_4_13 = 0.19758089021828817
# neg_flow_4_13 = 0.0


# time 6
# h_4 = 0.7681818181818182
# h_13 = 0.7681818181818182
# pos_flow_4_13 = 0.0
# neg_flow_4_13 = 0.0