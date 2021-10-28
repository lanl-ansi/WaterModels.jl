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
    set_flow_partitions_si!(subnetwork_mn, 10.0, 1.0e-4)
    wm = WaterModels.instantiate_model(subnetwork_mn, PWLRDWaterModel, build_owf_geo_part;)
    #relax_all_binary_variables!(wm)
    # JuMP.set_optimizer(wm.model, Xpress.Optimizer)
    # JuMP.set_optimizer_attributes(wm.model, "MAXTIME"=>300.0)
    # JuMP.set_optimizer_attributes(wm.model, "logfile" => "output.log")
    append!(wms, [wm])
end

wm = wms[2]
xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>30.0);

# pos_set_points = [0.13975120381145903, 0.11235517053879783, 0.1299985815974241, 0.15661160289348072, 0.112981175472513, 0.12036812776579112]
# neg_set_points = [8.852332289497124e-9, 6.631073181199603e-8, 5.185388917163471e-8, 3.450655071186831e-7, 5.622137396013869e-7, 3.894163388549295e-8]

pos_set_points = [0.1397566595881485, 0.11235861754645814, 0.12999822048682458, 0.15661495392588198, 0.11298094871356107, 0.12036937013952273]
neg_set_points = [8.852332289497124e-9, 6.631073181199603e-8, 5.185388917163471e-8, 3.450655071186831e-7, 5.622137396013869e-7, 3.894163388549295e-8]
for time in 1:(length(pos_set_points)) # 1 ok 3 ok 4 ok

    vref = wm.var[:it][:wm][:nw][time][:qp_pipe][2]
    JuMP.fix(vref, pos_set_points[time]; force = true)

    vref = wm.var[:it][:wm][:nw][time][:qn_pipe][2]
    JuMP.fix(vref, neg_set_points[time]; force = true)
end
result_LP = WaterModels.optimize_model!(wm, optimizer = xp)

# xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>30.0);
# result_LP = WaterModels.optimize_model!(wms[2], optimizer = xp)

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
