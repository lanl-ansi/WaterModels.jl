T_min_steam = 100
T_max_water = 100
T_ext = 25 #+273

## Objective function selection
obj_energy = 1
obj_flow_pressure = 0
@assert(obj_energy + obj_flow_pressure == 1)

# # Plant parameters
T_plant = 115
p_plant = 350000 #pascals #45 psi
# f_plant = 3*40*(pi/4)*(0.3048)^2; #2

# Load parameters
# ************
# *** Note: For convenience, for the time being, loads are edited in the input_x.jl file
# ************

# q_fixed_load = Dict()
# for e in LoadSet
#     # get!(q_fixed_load, e, loads_dict[e]["Q"])
# get!(q_fixed_load, e, 0)
#     get!(q_fixed_load, 14, (3.425e6)/1)
#     get!(q_fixed_load, 15, (3.425e6)/1)
#     get!(q_fixed_load, 16, (3.425e6)/1)
#     get!(q_fixed_load, 17, (3.425e6)/1)
#     get!(q_fixed_load, 18, 1*(3.425e6)/1)
# end


# Variable bounds and start points and tolerances
#temperature
T_lb = 80
T_ub = 150
T_start = 100

#pressure
# p_tol = 1.38e4 #pascal Tolerance

# Pressure lower bound for all the nodes in the network
# p_lb = 34000*ones(n_junctions) #Pa
p_lb = Dict()
for i in SetN
    get!(p_lb,i,34000)
end
#Pressure lower bound for plant outlet
for e in PlantSet
    get!(p_lb, N_to[e], 275000)
    # p_lb[N_to[e]] = 275000
end

# Warm starting values for pressures (for optimization)
# p_start = zeros(n_junctions)
p_start = Dict()

for i in SetNO
    get!(p_start, i, p_plant)
    # p_start[i] = p_plant
end
for i in SetNR
    get!(p_start, i, p_lb[i])
    # p_start[i] = p_lb[i]
end

#flow
phi_lb = 0
f_lb = 0
f_start = 3*40*(pi/4)*(0.3048)^2
# f_start = f_plant
