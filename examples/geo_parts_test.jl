using Pkg
Pkg.activate("/Users/dsigler/WaterModels.jl")

using Distributed

const WORKERS = 1 # Change to > 1 to use parallel
if nworkers() < WORKERS
    diff = (nprocs() == nworkers() ? WORKERS : WORKERS - nworkers())
    println("Adding $diff worker processes.")
    Distributed.addprocs(diff)
end

#Make sure these workers also have an environment with PH installed
@everywhere using Pkg
for w in workers()
    @spawnat(w, Pkg.activate("/Users/dsigler/WaterModels.jl"))
end

@everywhere using JuMP
#@everywhere using Ipopt
#@everywhere using TimerOutputs
@everywhere import Xpress
#@everywhere import Gurobi
@everywhere using Revise
@everywhere using WaterModels
#@everywhere using GLPK




network  = parse_file("/Users/dsigler/WaterModels.jl/examples/data/epanet/van_zyl.inp");
network_mn = WaterModels.make_multinetwork(network);

#ext = Dict{Symbol, Any}(:pipe_breakpoints => 5, :pump_breakpoints => 5)

#gx, node_map, w = build_graph_meta(network)

partition = [1, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 2, 2, 1, 1, 2]
br_pts    = [["1", "3", "7", "8", "10", "11", "13", "15", "16","17"],
                ["2", "4", "5", "6", "9", "12", "14"]]
# out = create_node_partition(network, 2)
# split_networks = split_multinetwork(network_mn,
# [["1", "2", "3", "4", "5","6","7"],
# ["7", "8","9","10","11","12","13"],
# ["13","14","15","16","17","18"],
# ["18","19","20","21","22","23","24"]]);

geosplit_networks = spatial_partition(network, br_pts)

##

wms = Array{AbstractWaterModel}([]) # Initialize empty array of WaterModels objects.
count = 1
# network_mn  = WaterModels.make_multinetwork(network)
# wm_ = WaterModels.instantiate_model(network_mn, PWLRDWaterModel, build_mn_owf; ext=ext)
for subnetwork in geosplit_networks
    #ext_ = Dict{Symbol, Any}(:pipe_breakpoints => 5, :pump_breakpoints => 5)
    subnetwork_mn = WaterModels.make_multinetwork(subnetwork)
    set_flow_partitions_si!(subnetwork_mn, 1.0, 1.0e-4)
    wm = WaterModels.instantiate_model(subnetwork_mn, PWLRDWaterModel, build_owf_geo_part;)
    #relax_all_binary_variables!(wm)
    # JuMP.set_optimizer(wm.model, Xpress.Optimizer)
    # JuMP.set_optimizer_attributes(wm.model, "MAXTIME"=>300.0)
    # JuMP.set_optimizer_attributes(wm.model, "logfile" => "output.log")
    append!(wms, [wm])
end

wm = wms[2]
xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>30.0);

part_1_pos_array = Any[0.22754379903928512, 0.23265535403121737, 0.23723108281298882, 0.1247561863541095]
#part_2_pos_array = Any[0.2017461603158504, 0.2119500071602848, 0.21115347700022907, 0.12433755975556912]
part_2_pos_array = Any[0.3700765647359651, 0.3700765647359655, 0.37007659606969445, 0.0]#
part_1_neg_array = Any[1.804665269884465e-8, -0.0, -0.0, 6.791945340594813e-8]
part_2_neg_array = Any[-0.0, -0.0, -0.0, -0.0]
pos_set_points = part_2_pos_array
neg_set_points = part_2_neg_array

# pos_flow_array = Any[0.2411775428041607, 0.2598190603862963, 0.19643361464907033, 0.0]
# neg_flow_array = Any[0.0, 0.0, 0.0, 0.0]
# h_4_array = Any[1.0247915558722451, 0.9894384162847942, 0.8635675455329749, 0.7681818181818182]
# h_13_array = Any[0.8793539521371608, 0.8225119671669396, 0.7697081113168751, 0.7681818181818182]
# pos_set_points = pos_flow_array
# neg_set_points = neg_flow_array



for time in 1:(length(pos_set_points)) # 1 ok 3 ok 4 ok

    vref = wm.var[:it][:wm][:nw][time][:qp_pipe][2]
    JuMP.fix(vref, pos_set_points[time]; force = true)

    vref = wm.var[:it][:wm][:nw][time][:qn_pipe][2]
    JuMP.fix(vref, neg_set_points[time]; force = true)

    # vref = wm.var[:it][:wm][:nw][time][:h][4]
    # JuMP.fix(vref, h_4_array[time]; force = true)

    # vref = wm.var[:it][:wm][:nw][time][:h][13]
    # JuMP.fix(vref, h_13_array[time]; force = true)

end
result_LP = WaterModels.optimize_model!(wm, optimizer = xp)

# xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>30.0);
# result_LP = WaterModels.optimize_model!(wms[2], optimizer = xp)

for time in  1:4
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

pipes = [2,10,11]
for time in  1:4
    println("")
    println("time ",time)
    for pipe in pipes
        @show pipe
        pos_flow = JuMP.value(wm.var[:it][:wm][:nw][time][:qp_pipe][pipe])
        neg_flow = JuMP.value(wm.var[:it][:wm][:nw][time][:qn_pipe][pipe])
        @show pos_flow
        @show neg_flow
    end


    println("")
end