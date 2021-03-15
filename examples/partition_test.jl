import Xpress
import JuMP
using Revise
using WaterModels
 
network = parse_file("examples/data/epanet/van_zyl.inp");
network_mn = WaterModels.make_multinetwork(network);
ext = Dict{Symbol, Any}(:pipe_breakpoints => 5, :pump_breakpoints => 5)
wm = WaterModels.instantiate_model(network_mn , PWLRDWaterModel, build_mn_owf; ext = ext)
relax_all_binary_variables!(wm)
xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>600.0);
result_LP = WaterModels.optimize_model!(wm, optimizer = xp)

#constraint_tank_volume_fixed(wm, n::Int, i::Int, V_0::Float64)

#split_networks = split_multinetwork(network_mn, [["1", "2"], ["2", "3", "4", "5"], ["6", "7", "8"]]);
#split_networks = split_multinetwork(network_mn, [["1", "2", "3", "4", "5","6", "7", "8"]]);
split_networks = split_multinetwork(network_mn, [["1", "2", "3", "4", "5","6"], ["6","7", "8","9","10","11","12"],
["12","13","14","15","16","17"],["17","18","19","20","21","22","23","24"]]);

# v_{t} = v_{t-1} + (flow_in_{t-1} - flow_out_{t-1})constants    

# v_{7} = v_{6} + (flow_in_{6} - flow_out_{6})constants  

# problem 1 has v_6

# somewhere there needs to be calcution to get v_7

# pass v_7 to problem 2 

# problem 2 has v_7

# v_7 is an input to problem 2 
 
wms = Array{AbstractWaterModel}([]) # Initialize empty array of WaterModels objects.
count = 1

for subnetwork in split_networks
    global count
    if count == 1
        ext = Dict{Symbol, Any}(:pipe_breakpoints => 5, :pump_breakpoints => 5)
        wm = WaterModels.instantiate_model(subnetwork, PWLRDWaterModel, build_mn_owf_part_start; ext = ext)
        append!(wms, [wm])
    else
        ext = Dict{Symbol, Any}(:pipe_breakpoints => 5, :pump_breakpoints => 5)
        wm = WaterModels.instantiate_model(subnetwork, PWLRDWaterModel, build_mn_owf_part; ext = ext)
        append!(wms, [wm])
    end
    count += 1
end

for prob in wms
    # time we care about at seems
    global_idx = sort([parse(Int, x) for x in collect( keys(prob.data["nw"]))])
    global_times = [first(global_idx),last(global_idx)]
    #local_times = [1,length(global_idx)]

    for global_time in global_times
        # get relevant global and local times
        #global_time = global_times[idx] 
        #local_time = local_times[idx]
        # loop over tanks
        for tank in keys(wm.data["nw"][string(global_time)]["tank"])
            # grab tank node
            node  = wm.data["nw"][string(global_time)]["tank"][tank]["node"]
            # grab LP tank head
            lp_var_ref = wm.var[:nw][global_time][:h][node]
            lp_value = JuMP.value(lp_var_ref)
            if tank =="2" 
                println("time",global_time)
                println("tank",tank)
                println(lp_var_ref,lp_value)
            end
            prob_var_ref = prob.var[:nw][global_time][:h][node]
            JuMP.fix(prob_var_ref,lp_value,force=true)
        end
    end
        

end


xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>300.0);

for part in wms

    println(" ")
    println(" ")
    println(" ")
    println(" ")
    println(" ")
    println(" ")
    println("starting new problem")
    println(" ")
    println(" ")
    println(" ")
    println(" ")
    println(" ")
    println(" ")


    xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>15.0);
 
    # JuMP.set_optimizer(part.model, Xpress.Optimizer)
    # JuMP.set_optimizer_attributes(part.model, "MAXTIME"=>300.0)
    # JuMP.optimize!(part.model)

    result = WaterModels.optimize_model!(part, optimizer = xp)
end

# JuMP.set_optimizer(wms[1].model, Xpress.Optimizer)
# JuMP.set_optimizer_attributes(wms[1].model, "MAXTIME"=>300.0)
# JuMP.optimize!(wms[1].model)



# ext = Dict{Symbol, Any}(:pipe_breakpoints => 5, :pump_breakpoints => 5)
# wm = WaterModels.instantiate_model(network_mn, PWLRDWaterModel, build_mn_owf; ext = ext)
# xp = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MAXTIME"=>60.0);
# result = WaterModels.optimize_model!(wm, optimizer = xp)

#.011570 53 seconds


# Goal: 
#     we want the cost savings if flexible loads where scheduled optimally to help the power grid run at a minimal cost

# LB is a lower bound on the best cost you could achieve by scheduling flexible loads. It might be weak that over estimates the possible savings

# LB solve FW with Binaries relaxed 

# drawback:
#     - loads you get are infeasible
#     - loads when paired with actualy UC decisions result in a high cost 
#         - wants everything on half of mingen 

# UB 1:

#     run UC with loads from the LB 

#     output: commitment decisions and Loads

#     drawback:
#     - loads you get are infeasible
#     - loads when paired with actualy UC decisions result in a high cost 
#         - wants everything on half of mingen 
#     pros: 
#     - fast 

# UB 2
#     some intial Loads # gives some binaries 
#     for i in range or until stopping criteria 
#         Solve UC with Loads
#         Fix binaries 
#         Run FW with binaries fixed until convergence
#         save FW_loads
#         Loads = FW_Loads

# output: commitment decisions and Loads

# drawbacks:
#         - could be a local minimum
#         - slower
# Pros:
#         you keep a feasible solution the whole time


# UB1-LB/min{LB,UB1} # this doesn't make sense unless we check that UB1 is feasible
# UB2-LB/min{LB,UB2}

# this might let you that your in within 2% of the best saving possible 

# example 

# Using UB2 we find a 1 million saving in the UC cost
# We compute LB and UB2-LB = 100K
# Then we can say that we have a solution that reschedules loads and give a 1 million savings. 
# Also we know that  the best possible savings are less than or equalt to 1,100,000  
# so our heurstic does a good job finds most of savings if not all. 

# power from gens at time t  == (y)newload_t - (z)oldload_t

# 1 = y + z