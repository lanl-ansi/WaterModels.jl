# Objective Functions

if(obj_energy == 1)
    # Minimize plant energy
    @objective(model, Min, sum(c_in[e]*f[e]*(100 - T_in[e]) + c_L*f[e] + c_out[e]*f[e]*(T_out[e] - 100) + p[N_fr[e]] for e in PlantSet))
elseif (obj_flow_pressure == 1)
    # Minimize plant flow and inlet pressure
    @objective(model, Min, sum(f[e] for e in PlantSet) + sum(p[N_fr[e]] for e in PlantSet))
end
