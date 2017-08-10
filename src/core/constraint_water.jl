##########################################################################################################
# The purpose of this file is to define commonly used and created constraints used in water flow models
##########################################################################################################


# Constraint that states a flow direction must be chosen
function constraint_flow_direction_choice{T}(wm::GenericWaterModel{T}, connection)
    i = connection["index"]

    yp = getvariable(wm.model, :yp)[i]
    yn = getvariable(wm.model, :yn)[i]

    c = @constraint(wm.model, yp + yn == 1)
    return Set([c]) # checked
end
#
# # Constraint that states a flow direction must be chosen
# function constraint_flow_direction_choice_ne{T}(wm::GenericWaterModel{T}, connection)
#     i = connection["index"]
#
#     yp = getvariable(wm.model, :yp_ne)[i]
#     yn = getvariable(wm.model, :yn_ne)[i]
#
#     c = @constraint(wm.model, yp + yn == 1)
#     return Set([c]) # checked
# end

function constraint_on_off_pipe_flow_bounds{T}(wm::GenericWaterModel{T}, pipe)
    pipe_idx = pipe["index"]
    i_junction_idx = pipe["f_junction"]
    j_junction_idx = pipe["t_junction"]

    yp = getvariable(wm.model, :yp)[pipe_idx]
    yn = getvariable(wm.model, :yn)[pipe_idx]
    f = getvariable(wm.model, :f)[pipe_idx]

    max_flow = wm.data["max_flow"]
    pd_max = pipe["pd_max"]
    pd_min = pipe["pd_min"]
    w = pipe["resistance"]

    c1 = @constraint(wm.model, -(1-yp)*min(max_flow, sqrt(w*max(pd_max, abs(pd_min)))) <= f)
    c2 = @constraint(wm.model, f <= (1-yn)*min(max_flow, sqrt(w*max(pd_max, abs(pd_min)))))

    return Set([c1, c2])    # checked (though one is reveresed of the other)
end

#######

# # constraints on flow across pipes where the directions are fixed
# function constraint_on_off_pipe_flow_bounds_fixed_direction{T}(wm::GenericWaterModel{T}, pipe)
#     pipe_idx = pipe["index"]
#     i_junction_idx = pipe["f_junction"]
#     j_junction_idx = pipe["t_junction"]
#
#     yp = pipe["yp"]
#     yn = pipe["yn"]
#     f = getvariable(wm.model, :f)[pipe_idx]
#
#     max_flow = wm.data["max_flow"]
#     pd_max = pipe["pd_max"]
#     pd_min = pipe["pd_min"]
#     w = pipe["resistance"]
#
#     c1 = @constraint(wm.model, -(1-yp)*min(max_flow, sqrt(w*max(pd_max, abs(pd_min)))) <= f)
#     c2 = @constraint(wm.model, f <= (1-yn)*min(max_flow, sqrt(w*max(pd_max, abs(pd_min)))))
#
#     return Set([c1, c2])    # checked (though one is reveresed of the other)
# end

# # constraints on flow across pipes
# function constraint_on_off_pipe_flow_direction_ne{T}(wm::GenericWaterModel{T}, pipe)
#     pipe_idx = pipe["index"]
#     i_junction_idx = pipe["f_junction"]
#     j_junction_idx = pipe["t_junction"]
#
#     yp = getvariable(wm.model, :yp_ne)[pipe_idx]
#     yn = getvariable(wm.model, :yn_ne)[pipe_idx]
#     f = getvariable(wm.model, :f_ne)[pipe_idx]
#
#     max_flow = wm.data["max_flow"]
#     pd_max = pipe["pd_max"]
#     pd_min = pipe["pd_min"]
#     w = pipe["resistance"]
#
#     c1 = @constraint(wm.model, -(1-yp)*min(max_flow, sqrt(w*max(pd_max, abs(pd_min)))) <= f)
#     c2 = @constraint(wm.model, f <= (1-yn)*min(max_flow, sqrt(w*max(pd_max, abs(pd_min)))))
#
#     return Set([c1, c2])    # checked (though one is reveresed of the other)
# end

# # constraints on flow across pipes when directions are fixed
# function constraint_on_off_pipe_flow_direction_ne_fixed_direction{T}(wm::GenericWaterModel{T}, pipe)
#     pipe_idx = pipe["index"]
#     i_junction_idx = pipe["f_junction"]
#     j_junction_idx = pipe["t_junction"]
#
#     yp = pipe["yp"]
#     yn = pipe["yn"]
#     f = getvariable(wm.model, :f_ne)[pipe_idx]
#
#     max_flow = wm.data["max_flow"]
#     pd_max = pipe["pd_max"]
#     pd_min = pipe["pd_min"]
#     w = pipe["resistance"]
#
#     c1 = @constraint(wm.model, -(1-yp)*min(max_flow, sqrt(w*max(pd_max, abs(pd_min)))) <= f)
#     c2 = @constraint(wm.model, f <= (1-yn)*min(max_flow, sqrt(w*max(pd_max, abs(pd_min)))))
#
#     return Set([c1, c2])    # checked (though one is reveresed of the other)
# end

#########
function constraint_on_off_potential_bounds{T}(wm::GenericWaterModel{T}, pipe)
    pipe_idx = pipe["index"]
    i_junction_idx = pipe["f_junction"]
    j_junction_idx = pipe["t_junction"]

    yp = getvariable(wm.model, :yp)[pipe_idx]
    yn = getvariable(wm.model, :yn)[pipe_idx]

    hi = getvariable(wm.model, :h_water)[i_junction_idx]
    hj = getvariable(wm.model, :h_water)[j_junction_idx]

    pd_min = pipe["pd_min"]
    pd_max = pipe["pd_max"]

    c1 = @constraint(wm.model, (1-yp) * pd_min <= hi - hj)
    c2 = @constraint(wm.model,hi - hj <= (1-yn)* pd_max)

    return Set([c1,c2])  # checked (though one is reversed of the other)
end

# # constraints on pressure drop across pipes
# function constraint_on_off_potential_bounds_fixed_direction{T}(wm::GenericWaterModel{T}, pipe)
#     pipe_idx = pipe["index"]
#     i_junction_idx = pipe["f_junction"]
#     j_junction_idx = pipe["t_junction"]
#
#     yp = pipe["yp"]
#     yn = pipe["yn"]
#
#     hi = getvariable(wm.model, :h_water)[i_junction_idx]
#     hj = getvariable(wm.model, :h_water)[j_junction_idx]
#
#     pd_min = pipe["pd_min"]
#     pd_max = pipe["pd_max"]
#
#     c1 = @constraint(wm.model, (1-yp) * pd_min <= hi - hj)
#     c2 = @constraint(wm.model, hi - hj <= (1-yn)* pd_max)
#
#     return Set([c1,c2])  # checked (though one is reversed of the other)
# end

# # constraints on pressure drop across pipes
# function constraint_on_off_pressure_drop_ne{T}(wm::GenericWaterModel{T}, pipe)
#     pipe_idx = pipe["index"]
#     i_junction_idx = pipe["f_junction"]
#     j_junction_idx = pipe["t_junction"]
#
#     yp = getvariable(wm.model, :yp_ne)[pipe_idx]
#     yn = getvariable(wm.model, :yn_ne)[pipe_idx]
#
#     hi = getvariable(wm.model, :h_water)[i_junction_idx]
#     hj = getvariable(wm.model, :h_water)[j_junction_idx]
#
#     pd_min = pipe["pd_min"]
#     pd_max = pipe["pd_max"]
#
#     c1 = @constraint(wm.model, (1-yp) * pd_min <= hi - hj)
#     c2 = @constraint(wm.model, hi - hj <= (1-yn)* pd_max)
#
#     return Set([c1,c2])  # checked (though one is reversed of the other)
# end


# # constraints on pressure drop across pipes when the direction is fixed
# function constraint_on_off_pressure_drop_ne_fixed_direction{T}(gm::GenericWaterModel{T}, pipe)
#     pipe_idx = pipe["index"]
#     i_junction_idx = pipe["f_junction"]
#     j_junction_idx = pipe["t_junction"]
#
#     yp = pipe["yp"]
#     yn = pipe["yn"]
#
#     hi = getvariable(wm.model, :h_water)[i_junction_idx]
#     hj = getvariable(wm.model, :h_water)[j_junction_idx]
#
#     pd_min = pipe["pd_min"]
#     pd_max = pipe["pd_max"]
#
#     c1 = @constraint(wm.model, (1-yp) * pd_min <= hi - hj)
#     c2 = @constraint(wm.model, hi - hj <= (1-yn)* pd_max)
#
#     return Set([c1,c2])  # checked (though one is reversed of the other)
# end




# standard flow balance equation where demand and production is fixed
function constraint_junction_flow_balance{T}(wm::GenericWaterModel{T}, junction)
    i = junction["index"]
    junction_branches = wm.set.junction_connections[i]

    f_branches = collect(keys(filter( (a, connection) -> connection["f_junction"] == i, wm.set.connections)))
    t_branches = collect(keys(filter( (a, connection) -> connection["t_junction"] == i, wm.set.connections)))

    h = getvariable(wm.model, :h_water)
    f = getvariable(wm.model, :f)

    c = @constraint(wm.model, junction["qgfirm"] - junction["qlfirm"] == sum(f[a] for a in f_branches) - sum(f[a] for a in t_branches) )

    return Set([c])
end

# # standard flow balance equation where demand and production is fixed
# function constraint_junction_flow_balance_ne{T}(wm::GenericWaterModel{T}, junction)
#     i = junction["index"]
#     junction_branches = wm.set.junction_connections[i]
#
#     f_branches = collect(keys(filter( (a, connection) -> connection["f_junction"] == i, wm.set.connections)))
#     t_branches = collect(keys(filter( (a, connection) -> connection["t_junction"] == i, wm.set.connections)))
#
#     f_branches_ne = collect(keys(filter( (a, connection) -> connection["f_junction"] == i, wm.set.new_connections)))
#     t_branches_ne = collect(keys(filter( (a, connection) -> connection["t_junction"] == i, wm.set.new_connections)))
#
#     h = getvariable(wm.model, :h_water)
#     f = getvariable(wm.model, :f)
#     f_ne = getvariable(wm.model, :f_ne)
#     c = @constraint(wm.model, junction["qgfirm"] - junction["qlfirm"] == sum(f[a] for a in f_branches) - sum(f[a] for a in t_branches) + sum(f_ne[a] for a in f_branches_ne) - sum(f_ne[a] for a in t_branches_ne) )
#
#     return Set([c])
# end


# # standard flow balance equation where demand and production is fixed
# function constraint_junction_flow_balance_ls{T}(wm::GenericWaterModel{T}, junction)
#     i = junction["index"]
#     junction_branches = wm.set.junction_connections[i]
#
#     f_branches = collect(keys(filter( (a, connection) -> connection["f_junction"] == i, wm.set.connections)))
#     t_branches = collect(keys(filter( (a, connection) -> connection["t_junction"] == i, wm.set.connections)))
#
#     p = getvariable(wm.model, :p_gas)
#     f = getvariable(wm.model, :f)
#     ql = 0
#     qg = 0
#     if junction["qlmin"] != junction["qlmax"]
#         ql = getvariable(wm.model, :ql_gas)[i]
#     end
#     if junction["qgmin"] != junction["qgmax"]
#         qg = getvariable(wm.model, :qg_gas)[i]
#     end
#     ql_firm = junction["qlfirm"]
#     qg_firm = junction["qgfirm"]
#
#     c = @constraint(wm.model, qg_firm - ql_firm + qg - ql == sum(f[a] for a in f_branches) - sum(f[a] for a in t_branches) )
#
#     return Set([c])
# end

# standard flow balance equation where demand and production is fixed
# function constraint_junction_flow_balance_ne_ls{T}(wm::GenericWaterModel{T}, junction)
#     i = junction["index"]
#     junction_branches = wm.set.junction_connections[i]
#
#     f_branches = collect(keys(filter( (a, connection) -> connection["f_junction"] == i, wm.set.connections)))
#     t_branches = collect(keys(filter( (a, connection) -> connection["t_junction"] == i, wm.set.connections)))
#
#     f_branches_ne = collect(keys(filter( (a, connection) -> connection["f_junction"] == i, wm.set.new_connections)))
#     t_branches_ne = collect(keys(filter( (a, connection) -> connection["t_junction"] == i, wm.set.new_connections)))
#
#     p = getvariable(wm.model, :p_gas)
#     f = getvariable(wm.model, :f)
#     f_ne = getvariable(wm.model, :f_ne)
#
#     ql = 0
#     qg = 0
#     if junction["qlmin"] != junction["qlmax"]
#         ql = getvariable(wm.model, :ql_gas)[i]
#     end
#     if junction["qgmin"] != junction["qgmax"]
#         qg = getvariable(wm.model, :qg_gas)[i]
#     end
#
#     ql_firm = junction["qlfirm"]
#     qg_firm = junction["qgfirm"]
#
#     c = @constraint(wm.model, qg_firm - ql_firm + qg - ql == sum(f[a] for a in f_branches) - sum(f[a] for a in t_branches) + sum(f_ne[a] for a in f_branches_ne) - sum(f_ne[a] for a in t_branches_ne) )
#
#     return Set([c])
# end



# **************************************************************************************************************************


function constraint_head_flow_coupling{T}(wm::GenericWaterModel{T},pipe)
  pipe_idx = pipe["index"]
  i_junction_idx = pipe["f_junction"]
  j_junction_idx = pipe["t_junction"]

  yp = getvariable(wm.model, :yp)[pipe_idx]
  yn = getvariable(wm.model, :yn)[pipe_idx]
  f = getvariable(wm.model, :f)[pipe_idx]

  hi = getvariable(wm.model, :h_water)[i_junction_idx]
  hj = getvariable(wm.model, :h_water)[j_junction_idx]

  c1 = constraint(wm.model,)



#********************************************************************************************************************************


# constraints on flow across short pipes
function constraint_on_off_short_pipe_flow_direction{T}(wm::GenericWaterModel{T}, pipe)
    pipe_idx = pipe["index"]
    i_junction_idx = pipe["f_junction"]
    j_junction_idx = pipe["t_junction"]

    yp = getvariable(wm.model, :yp)[pipe_idx]
    yn = getvariable(wm.model, :yn)[pipe_idx]
    f = getvariable(wm.model, :f)[pipe_idx]
    max_flow = wm.data["max_flow"]

    c1 = @constraint(wm.model, -max_flow*(1-yp) <= f)
    c2 = @constraint(wm.model, f <= max_flow*(1-yn))

    return Set([c1, c2])
end

# # constraints on flow across short pipes when the directions are constants
# function constraint_on_off_short_pipe_flow_direction_fixed_direction{T}(wm::GenericWaterModel{T}, pipe)
#     pipe_idx = pipe["index"]
#     i_junction_idx = pipe["f_junction"]
#     j_junction_idx = pipe["t_junction"]
#
#     yp = pipe["yp"]
#     yn = pipe["yn"]
#     f = getvariable(wm.model, :f)[pipe_idx]
#     max_flow = wm.data["max_flow"]
#
#     c1 = @constraint(wm.model, -max_flow*(1-yp) <= f)
#     c2 = @constraint(wm.model, f <= max_flow*(1-yn))
#
#     return Set([c1, c2])
# end

# constraints on pressure drop across pipes
function constraint_short_pipe_pressure_drop{T}(wm::GenericWaterModel{T}, pipe)
    pipe_idx = pipe["index"]
    i_junction_idx = pipe["f_junction"]
    j_junction_idx = pipe["t_junction"]

    hi = getvariable(wm.model, :h_water)[i_junction_idx]
    hj = getvariable(wm.model, :h_water)[j_junction_idx]

    c = @constraint(wm.model,  pi == pj)
    return Set([c])
end

# constraints on flow across valves
function constraint_on_off_valve_flow_direction{T}(wm::GenericWaterModel{T}, valve)
    valve_idx = valve["index"]
    i_junction_idx = valve["f_junction"]
    j_junction_idx = valve["t_junction"]

    yp = getvariable(wm.model, :yp)[valve_idx]
    yn = getvariable(wm.model, :yn)[valve_idx]
    f = getvariable(wm.model, :f)[valve_idx]
    v = getvariable(wm.model, :valve)[valve_idx]

    max_flow = wm.data["max_flow"]

    c1 = @constraint(wm.model, -max_flow*(1-yp) <= f <= max_flow*(1-yn))
    c2 = @constraint(wm.model, -max_flow*v <= f <= max_flow*v)

    return Set([c1, c2])
end

# # constraints on flow across valves when directions are constants
# function constraint_on_off_valve_flow_direction_fixed_direction{T}(wm::GenericWaterModel{T}, valve)
#     valve_idx = valve["index"]
#     i_junction_idx = valve["f_junction"]
#     j_junction_idx = valve["t_junction"]
#
#     yp = valve["yp"]
#     yn = valve["yn"]
#     f = getvariable(wm.model, :f)[valve_idx]
#     v = getvariable(wm.model, :valve)[valve_idx]
#
#     max_flow = wm.data["max_flow"]
#
#     c1 = @constraint(wm.model, -max_flow*(1-yp) <= f <= max_flow*(1-yn))
#     c2 = @constraint(wm.model, -max_flow*v <= f <= max_flow*v)
#
#     return Set([c1, c2])
# end

# constraints on pressure drop across valves
function constraint_on_off_valve_pressure_drop{T}(wm::GenericWaterModel{T}, valve)
    valve_idx = valve["index"]
    i_junction_idx = valve["f_junction"]
    j_junction_idx = valve["t_junction"]

    i = gm.set.junctions[i_junction_idx]
    j = gm.set.junctions[j_junction_idx]

    hi = getvariable(wm.model, :h_water)[i_junction_idx]
    hj = getvariable(wm.model, :h_water)[j_junction_idx]

    v = getvariable(wm.model, :valve)[valve_idx]

    c1 = @constraint(wm.model,  hj - ((1-v)*j["pmax"]^2) <= hi)
    c2 = @constraint(wm.model,  hi <= hj + ((1-v)*i["pmax"]^2))

    return Set([c1, c2])
end

# constraints on flow across control valves
function constraint_on_off_control_valve_flow_direction{T}(wm::GenericWaterModel{T}, valve)
    valve_idx = valve["index"]
    i_junction_idx = valve["f_junction"]
    j_junction_idx = valve["t_junction"]

    yp = getvariable(wm.model, :yp)[valve_idx]
    yn = getvariable(wm.model, :yn)[valve_idx]
    f = getvariable(wm.model, :f)[valve_idx]
    v = getvariable(wm.model, :valve)[valve_idx]
    max_flow = wm.data["max_flow"]

    c1 = @constraint(wm.model, -max_flow*(1-yp) <= f)
    c2 = @constraint(wm.model, f <= max_flow*(1-yn))
    c3 = @constraint(wm.model, -max_flow*v <= f )
    c4 = @constraint(wm.model, f <= max_flow*v)

    return Set([c1, c2, c3, c4])
end

# constraints on flow across control valves when directions are constants
function constraint_on_off_control_valve_flow_direction_fixed_direction{T}(wm::GenericWaterModel{T}, valve)
    valve_idx = valve["index"]
    i_junction_idx = valve["f_junction"]
    j_junction_idx = valve["t_junction"]

    yp = valve["yp"]
    yn = valve["yn"]
    f = getvariable(wm.model, :f)[valve_idx]
    v = getvariable(wm.model, :valve)[valve_idx]
    max_flow = wm.data["max_flow"]

    c1 = @constraint(wm.model, -max_flow*(1-yp) <= f)
    c2 = @constraint(wm.model, f <= max_flow*(1-yn))
    c3 = @constraint(wm.model, -max_flow*v <= f )
    c4 = @constraint(wm.model, f <= max_flow*v)

    return Set([c1, c2, c3, c4])
end


# constraints on pressure drop across control valves
function constraint_on_off_control_valve_pressure_drop{T}(wm::GenericWaterModel{T}, valve)
    valve_idx = valve["index"]
    i_junction_idx = valve["f_junction"]
    j_junction_idx = valve["t_junction"]

    i = gm.set.junctions[i_junction_idx]
    j = gm.set.junctions[j_junction_idx]

    hi = getvariable(wm.model, :h_water)[i_junction_idx]
    hj = getvariable(wm.model, :h_water)[j_junction_idx]
    yp = getvariable(wm.model, :yp)[valve_idx]
    yn = getvariable(wm.model, :yn)[valve_idx]
    v = getvariable(wm.model, :valve)[valve_idx]

    max_ratio = valve["c_ratio_max"]
    min_ratio = valve["c_ratio_min"]

    c1 = @constraint(wm.model,  hj - ((1-v)*j["hmax"]) <= hi)
    c2 = @constraint(wm.model,  hi <= hj + ((1-v)*i["hmax"]))

    return Set([c1, c2])
end

# # constraints on pressure drop across control valves when directions are constants
# function constraint_on_off_control_valve_pressure_drop_fixed_direction{T}(gm::GenericGasModel{T}, valve)
#     valve_idx = valve["index"]
#     i_junction_idx = valve["f_junction"]
#     j_junction_idx = valve["t_junction"]
#
#     i = gm.set.junctions[i_junction_idx]
#     j = gm.set.junctions[j_junction_idx]
#
#     pi = getvariable(gm.model, :p_gas)[i_junction_idx]
#     pj = getvariable(gm.model, :p_gas)[j_junction_idx]
#     yp = valve["yp"]
#     yn = valve["yn"]
#     v = getvariable(gm.model, :valve)[valve_idx]
#
#     max_ratio = valve["c_ratio_max"]
#     min_ratio = valve["c_ratio_min"]
#
#     c1 = @constraint(gm.model,  pj - (max_ratio*pi) <= (2-yp-v)*j["pmax"]^2)
#     c2 = @constraint(gm.model,  (min_ratio*pi) - pj <= (2-yp-v)*(min_ratio*i["pmin"]^2) )
#     c3 = @constraint(gm.model,  pi - (max_ratio*pj) <= (2-yn-v)*i["pmax"]^2)
#     c4 = @constraint(gm.model,  (min_ratio*pj) - pi <= (2-yn-v)*(min_ratio*j["pmin"]^2))
#
#     return Set([c1, c2, c3, c4])
# end


# Make sure there is at least one direction set to take flow away from a junction (typically used on source nodes)
function constraint_source_flow{T}(wm::GenericWaterModel{T}, junction)
    f_branches = collect(keys(filter( (a,connection) -> connection["f_junction"] == junction["index"], wm.set.connections)))
    t_branches = collect(keys(filter( (a,connection) -> connection["t_junction"] == junction["index"], wm.set.connections)))

    yp = getvariable(wm.model, :yp)
    yn = getvariable(wm.model, :yn)

    c = @constraint(wm.model, sum(yp[a] for a in f_branches) + sum(yn[a] for a in t_branches) >= 1)
    return Set([c])
end

# Make sure there is at least one direction set to take flow away from a junction (typically used on source nodes)
# function constraint_source_flow_ne{T}(gm::GenericGasModel{T}, junction)
#     f_branches = collect(keys(filter( (a,connection) -> connection["f_junction"] == junction["index"], gm.set.connections)))
#     t_branches = collect(keys(filter( (a,connection) -> connection["t_junction"] == junction["index"], gm.set.connections)))
#
#     f_branches_ne = collect(keys(filter( (a,connection) -> connection["f_junction"] == junction["index"], gm.set.new_connections)))
#     t_branches_ne = collect(keys(filter( (a,connection) -> connection["t_junction"] == junction["index"], gm.set.new_connections)))
#
#     yp = getvariable(gm.model, :yp)
#     yn = getvariable(gm.model, :yn)
#
#     yp_ne = getvariable(gm.model, :yp_ne)
#     yn_ne = getvariable(gm.model, :yn_ne)
#
#     c = @constraint(gm.model, sum(yp[a] for a in f_branches) + sum(yn[a] for a in t_branches) + sum(yp_ne[a] for a in f_branches_ne) + sum(yn_ne[a] for a in t_branches_ne) >= 1)
#     return Set([c])
# end

# Make sure there is at least one direction set to take flow to a junction (typically used on sink nodes)
function constraint_sink_flow{T}(wm::GenericWaterModel{T}, junction)
    f_branches = collect(keys(filter( (a,connection) -> connection["f_junction"] == junction["index"], wm.set.connections)))
    t_branches = collect(keys(filter( (a,connection) -> connection["t_junction"] == junction["index"], wm.set.connections)))

    yp = getvariable(wm.model, :yp)
    yn = getvariable(wm.model, :yn)

    c = @constraint(wm.model, sum(yn[a] for a in f_branches) + sum(yp[a] for a in t_branches) >= 1)
    return Set([c])
end

# # Make sure there is at least one direction set to take flow to a junction (typically used on sink nodes)
# function constraint_sink_flow_ne{T}(gm::GenericGasModel{T}, junction)
#     f_branches = collect(keys(filter( (a,connection) -> connection["f_junction"] == junction["index"], gm.set.connections)))
#     t_branches = collect(keys(filter( (a,connection) -> connection["t_junction"] == junction["index"], gm.set.connections)))
#     f_branches_ne = collect(keys(filter( (a,connection) -> connection["f_junction"] == junction["index"], gm.set.new_connections)))
#     t_branches_ne = collect(keys(filter( (a,connection) -> connection["t_junction"] == junction["index"], gm.set.new_connections)))
#
#     yp = getvariable(gm.model, :yp)
#     yn = getvariable(gm.model, :yn)
#     yp_ne = getvariable(gm.model, :yp_ne)
#     yn_ne = getvariable(gm.model, :yn_ne)
#
#     c = @constraint(gm.model, sum(yn[a] for a in f_branches) + sum(yp[a] for a in t_branches) + sum(yn_ne[a] for a in f_branches_ne) + sum(yp_ne[a] for a in t_branches_ne) >= 1)
#     return Set([c])
# end

# This constraint is intended to ensure that flow is on direction through a node with degree 2 and no production or consumption
function constraint_conserve_flow{T}(gm::GenericGasModel{T}, junction)
    idx = junction["index"]
    first = nothing
    last = nothing

    for i in gm.set.junction_connections[idx]
        connection = gm.set.connections[i]
        if connection["f_junction"] == idx
            other = connection["t_junction"]
        else
            other = connection["f_junction"]
        end

        if first == nothing
            first = other
        elseif first != other
            last = other
        end
    end

    yp_first = filter(i -> gm.set.connections[i]["f_junction"] == first, gm.set.junction_connections[idx])
    yn_first = filter(i -> gm.set.connections[i]["t_junction"] == first, gm.set.junction_connections[idx])
    yp_last  = filter(i -> gm.set.connections[i]["t_junction"] == last,  gm.set.junction_connections[idx])
    yn_last  = filter(i -> gm.set.connections[i]["f_junction"] == last,  gm.set.junction_connections[idx])

    yp = getvariable(gm.model, :yp)
    yn = getvariable(gm.model, :yn)

    c =  Set([])
    if length(yn_first) > 0 && length(yp_last) > 0
        for i1 in yn_first
            for i2 in yp_last
              c1 = @constraint(gm.model, yn[i1]  == yp[i2])
              c = Set([c;c1])
            end
        end
    end


   if length(yn_first) > 0 && length(yn_last) > 0
        for i1 in yn_first
            for i2 in yn_last
              c1 = @constraint(gm.model, yn[i1] == yn[i2])
              c = Set([c;c1])
            end
        end
    end

    if length(yp_first) > 0 && length(yp_last) > 0
        for i1 in yp_first
            for i2 in yp_last
              c1 = @constraint(gm.model, yp[i1]  == yp[i2])
              c = Set([c;c1])
            end
        end
    end

    if length(yp_first) > 0 && length(yn_last) > 0
        for i1 in yp_first
            for i2 in yn_last
              c1 = @constraint(gm.model, yp[i1] == yn[i2])
              c = Set([c;c1])
            end
        end
    end

    return c
end


# # This constraint is intended to ensure that flow is on direction through a node with degree 2 and no production or consumption
# function constraint_conserve_flow_ne{T}(gm::GenericGasModel{T}, junction)
#     idx = junction["index"]
#     first = nothing
#     last = nothing
#
#     for i in gm.set.junction_connections[idx]
#         connection = gm.set.connections[i]
#         if connection["f_junction"] == idx
#             other = connection["t_junction"]
#         else
#             other = connection["f_junction"]
#         end
#
#         if first == nothing
#             first = other
#         elseif first != other
#             last = other
#         end
#     end
#
#     for i in gm.set.junction_new_connections[idx]
#         connection = gm.set.new_connections[i]
#         if connection["f_junction"] == idx
#             other = connection["t_junction"]
#         else
#             other = connection["f_junction"]
#         end
#
#         if first == nothing
#             first = other
#         elseif first != other
#             last = other
#         end
#     end
#
#
#     yp_first = [filter(i -> gm.set.connections[i]["f_junction"] == first, gm.set.junction_connections[idx]); filter(i -> gm.set.new_connections[i]["f_junction"] == first, gm.set.junction_new_connections[idx])]
#     yn_first = [filter(i -> gm.set.connections[i]["t_junction"] == first, gm.set.junction_connections[idx]); filter(i -> gm.set.new_connections[i]["t_junction"] == first, gm.set.junction_new_connections[idx])]
#     yp_last  = [filter(i -> gm.set.connections[i]["t_junction"] == last,  gm.set.junction_connections[idx]); filter(i -> gm.set.new_connections[i]["t_junction"] == last,  gm.set.junction_new_connections[idx])]
#     yn_last  = [filter(i -> gm.set.connections[i]["f_junction"] == last,  gm.set.junction_connections[idx]); filter(i -> gm.set.new_connections[i]["f_junction"] == last,  gm.set.junction_new_connections[idx])]
#
#     yp = getvariable(gm.model, :yp)
#     yn = getvariable(gm.model, :yn)
#     yp_ne = getvariable(gm.model, :yp_ne)
#     yn_ne = getvariable(gm.model, :yn_ne)
#
#     c =  Set([])
#     if length(yn_first) > 0 && length(yp_last) > 0
#         for i1 in yn_first
#             for i2 in yp_last
#                 yn1 = haskey(gm.set.connections,i1) ? yn[i1] : yn_ne[i1]
#                 yp2 = haskey(gm.set.connections,i2) ? yp[i2] : yp_ne[i2]
#                 c1 = @constraint(gm.model, yn1  == yp2)
#                 c = Set([c;c1])
#             end
#         end
#     end
#
#     if length(yn_first) > 0 && length(yn_last) > 0
#         for i1 in yn_first
#             for i2 in yn_last
#                 yn1 = haskey(gm.set.connections,i1) ? yn[i1] : yn_ne[i1]
#                 yn2 = haskey(gm.set.connections,i2) ? yn[i2] : yn_ne[i2]
#                 c1 = @constraint(gm.model, yn1 == yn2)
#                 c = Set([c;c1])
#             end
#         end
#     end
#
#     if length(yp_first) > 0 && length(yp_last) > 0
#         for i1 in yp_first
#             for i2 in yp_last
#                 yp1 = haskey(gm.set.connections,i1) ? yp[i1] : yp_ne[i1]
#                 yp2 = haskey(gm.set.connections,i2) ? yp[i2] : yp_ne[i2]
#                 c1 = @constraint(gm.model, yp1  == yp2)
#                 c = Set([c;c1])
#             end
#         end
#     end
#
#     if length(yp_first) > 0 && length(yn_last) > 0
#         for i1 in yp_first
#             for i2 in yn_last
#                 yp1 = haskey(gm.set.connections,i1) ? yp[i1] : yp_ne[i1]
#                 yn2 = haskey(gm.set.connections,i2) ? yn[i2] : yn_ne[i2]
#                 c1 = @constraint(gm.model, yp1 == yn2)
#                 c = Set([c;c1])
#             end
#         end
#     end
#
#     return c
# end

# # ensures that parallel lines have flow in the same direction
# function constraint_parallel_flow{T}(gm::GenericGasModel{T}, connection)
#     i = min(connection["f_junction"], connection["t_junction"])
#     j = max(connection["f_junction"], connection["t_junction"])
#     idx = connection["index"]
#
#     f_connections = filter(i -> gm.set.connections[i]["f_junction"] == connection["f_junction"], gm.set.parallel_connections[(i,j)])
#     t_connections = filter(i -> gm.set.connections[i]["f_junction"] != connection["f_junction"], gm.set.parallel_connections[(i,j)])
#
#     yp = getvariable(gm.model, :yp)
#     yn = getvariable(gm.model, :yn)
#
#     c = @constraint(gm.model, sum(yp[i] for i in f_connections) + sum(yn[i] for i in t_connections) == yp[idx] * length(gm.set.parallel_connections[(i,j)]))
#     return Set([c])
# end


# # ensures that parallel lines have flow in the same direction
# function constraint_parallel_flow_ne{T}(gm::GenericGasModel{T}, connection)
#     i = min(connection["f_junction"], connection["t_junction"])
#     j = max(connection["f_junction"], connection["t_junction"])
#     idx = connection["index"]
#
#     f_connections = filter(i -> gm.set.connections[i]["f_junction"] == connection["f_junction"], intersect(gm.set.all_parallel_connections[(i,j)], gm.set.parallel_connections[(i,j)]))
#     t_connections = filter(i -> gm.set.connections[i]["f_junction"] != connection["f_junction"], intersect(gm.set.all_parallel_connections[(i,j)], gm.set.parallel_connections[(i,j)]))
#     f_connections_ne = filter(i -> gm.set.new_connections[i]["f_junction"] == connection["f_junction"], setdiff(gm.set.all_parallel_connections[(i,j)], gm.set.parallel_connections[(i,j)]))
#     t_connections_ne = filter(i -> gm.set.new_connections[i]["f_junction"] != connection["f_junction"], setdiff(gm.set.all_parallel_connections[(i,j)], gm.set.parallel_connections[(i,j)]))
#
#     yp = getvariable(gm.model, :yp)
#     yn = getvariable(gm.model, :yn)
#     yp_ne = getvariable(gm.model, :yp_ne)
#     yn_ne = getvariable(gm.model, :yn_ne)
#     yp_i = haskey(gm.set.connections, idx) ? yp[idx] : yp_ne[idx]
#
#     c = @constraint(gm.model, sum(yp[i] for i in f_connections) + sum(yn[i] for i in t_connections) + sum(yp_ne[i] for i in f_connections_ne) + sum(yn_ne[i] for i in t_connections_ne) == yp_i * length(gm.set.all_parallel_connections[(i,j)]))
#     return Set([c])
# end
