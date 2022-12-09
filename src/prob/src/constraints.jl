#
# # Nodal mass flow balance
@constraint(model, mass_flow_balance[i in setdiff(SetN)], sum(phi[ei]*A[ei] for ei in P_in[i]) +
sum(f[ei] for ei in L_in[i]) +
sum(f[ei] for ei in PL_in[i])
 == sum(phi[eo]*A[eo] for eo in P_out[i]) +
 sum(f[eo] for eo in L_out[i])+
 sum(f[eo] for eo in PL_out[i])
 )
#
# #Nodel thermal energy balance
@constraint(model, thermal_energy_balance[i in SetN],
sum(phi[e]*A[e]*c[e]*T_out[e] for e in P_in[i]) +
sum(f[e]*c_out[e]*T_out[e] for e in L_in[i]) +
sum(f[e]*c_out[e]*T_out[e] for e in PL_in[i])
==
sum(phi[e]*A[e]*c[e]*T_in[e] for e in P_out[i]) +
sum(f[e]*c_in[e]*T_in[e] for e in L_out[i]) +
sum(f[e]*c_in[e]*T_in[e] for e in PL_out[i])
)
# # Mixing at nodes
@constraint(model, mixing_at_nodes[i in SetN, e in setdiff(E_out[i])], T_in[e] == T[i])

# # #Thermal loss in outgoing/steam pipes
@NLconstraint(model, temperature_loss_steam[e in SetPO], T_out[e] == T_ext + (T_in[e] - T_ext)*exp(-L[e]/(c[e]*gamma[e]*phi[e]*A[e])))

# #Thermal loss in incoming/water pipes
# @NLconstraint(model, temperature_loss_water[e in SetPR], T_out[e] == T_ext + (T_in[e] - T_ext)*exp(-L[e]/(c[e]*gamma[e]*phi[e]*A[e])))
@NLconstraint(model, temperature_loss_water[e in SetPR], T_out[e] == T_ext + (T_in[e] - T_ext)*exp(-gamma[e]*L[e]*rho[e]/phi[e]))
#
# # #Pressure loss steam
@NLconstraint(model, pressure_loss_steam[e in SetPO], p[N_to[e]]^2 - p[N_fr[e]]^2 == (-lambda[e]*R_s/d[e])*phi[e]*abs(phi[e])*(T_ext*L[e] + c[e]*phi[e]*A[e]*gamma[e]*(T_in[e] - T_out[e])))

# # #Pressure loss water
@NLconstraint(model, pressure_loss_water[e in SetPR], p[N_fr[e]] - p[N_to[e]] == phi[e]*abs(phi[e])*lambda[e]*L[e]/(2*d[e]*rho[e]))

# #Plant
# @constraint(model, plant_temperature[e in PlantSet], T_out[e] == T_plant)
# @constraint(model, plant_pressure[e in PlantSet], p[N_to[e]] == p_plant)
# @constraint(model, plant_flux[e in PlantSet], f[e] == f_plant)

# #
# #Load
@NLconstraint(model, load_temperature_drop[e in LoadSet], c_in[e]*f[e]*(T_in[e]-100) + c_L*f[e] + c_out[e]*f[e]*(100-T_out[e]) >= Q_load[e])
@constraint(model, load_pressure_with_pipe_after_condensor_1[e in LoadSet], p[N_fr[e]] >= p[N_to[e]])
@constraint(model, steam_entering_load[e in LoadSet], T_in[e]>=T_min_steam)
@constraint(model, water_returning_from_load[e in LoadSet], T_out[e]<=T_max_water)
