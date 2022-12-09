c_in = Dict(); c_out = Dict()
for e in PlantSet
    get!(c_in,e,c_v_w)
    get!(c_out,e,c_p_s)
end

# # Load
#using a fixed value for all loads
Q_load = Dict()
for e in LoadSet
    # get!(Q_load, e, q_fixed_load[e])
    get!(Q_load, e, loads_dict[e]["Q"])
    get!(c_in,e,c_p_s)
    get!(c_out,e,c_v_w)
end

# # Pipe
c = Dict(); gamma = Dict(); lambda = Dict();
L = Dict(); d = Dict(); A = Dict(); rho = Dict();

# # Steam pipe
for e in SetPO
    get!(c,e,c_p_s)
    get!(gamma,e,pipes_dict[e]["thermal_loss_coefficient"])
    get!(lambda,e,pipes_dict[e]["friction_factor"])

    get!(L,e, pipes_dict[e]["length"])
    get!(d,e,pipes_dict[e]["diameter"])
    area = (pi/4)*d[e]^2;
    get!(A,e,area)
    get!(rho, e, rho_s)
end

# # Water pipe
for e in SetPR
    get!(c,e,c_v_w)
    get!(gamma,e,pipes_dict[e]["thermal_loss_coefficient"])
    get!(lambda,e,pipes_dict[e]["friction_factor"])

    get!(L, e, pipes_dict[e]["length"])
    get!(d,e,pipes_dict[e]["diameter"])
    area = (pi/4)*d[e]^2;
    get!(A,e,area)
    get!(rho, e, rho_w)
end
