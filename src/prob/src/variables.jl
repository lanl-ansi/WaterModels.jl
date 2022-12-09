# Variables
#junctions
@variable(model, T[i in SetN], lower_bound = T_lb, upper_bound = T_ub)
@variable(model, p[i in SetN], lower_bound = p_lb[i], start = p_start[i])

@variable(model, T_in[e in SetE],lower_bound = T_lb, upper_bound = T_ub, start = T_start)
@variable(model, T_out[e in SetE],lower_bound = T_lb, upper_bound = T_ub, start = T_start)

@variable(model, phi[e in SetP], lower_bound = phi_lb, start = f_start/A[e])
@variable(model, f[e in union(LoadSet, PlantSet)], lower_bound = f_lb,start = f_start)
