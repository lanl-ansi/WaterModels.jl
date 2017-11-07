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
#=
#Constraint that states a flow direction must be chosen for new edges
function constraint_flow_direction_choice_ne{T}(wm::GenericWaterModel{T}, i)
    yp = gm.var[:yp_ne][i]
    yn = gm.var[:yn_ne][i]

    c = @constraint(gm.model, yp + yn == 1)
    return Set([c]) # checked
end
=#
#constraints on pressure drop across pipes
function constraint_on_off_head_bounds{T}(wm::GenericWaterModel{T}, pipe)
    pipe_idx = pipe["index"]
    i_junction_idx = pipe["f_junction"]
    j_junction_idx = pipe["t_junction"]

    yp = getvariable(wm.model, :yp)[pipe_idx]
    yn = getvariable(wm.model, :yn)[pipe_idx]

    hi = getvariable(wm.model, :h_water)[i_junction_idx]
    hj = getvariable(wm.model, :h_water)[j_junction_idx]

    hd_min = pipe["hd_min"]
    hd_max = pipe["hd_max"]

    c1 = @constraint(wm.model, (1-yp) * hd_min <= hi - hj)
    c2 = @constraint(wm.model,hi - hj <= (1-yn)* hd_max)

    return Set([c1,c2])  # checked (though one is reversed of the other)
end

# #constraints on pressure drop across pipes
# function constraint_on_off_head_bounds_fixed_direction{T}(wm::GenericWaterModel{T}, pipe)
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
#     hd_min = pipe["hd_min"]
#     hd_max = pipe["hd_max"]
#
#     c1 = @constraint(wm.model, (1-yp) * hd_min <= hi - hj)
#     c2 = @constraint(wm.model,hi - hj <= (1-yn)* hd_max)
#
#     return Set([c1,c2])  # checked (though one is reversed of the other)
# end

#Bounds on pipe flow (by max flow overall)
function constraint_on_off_pipe_flow_bounds_max_flow{T}(wm::GenericWaterModel{T}, pipe)
    pipe_idx = pipe["index"]
    i_junction_idx = pipe["f_junction"]
    j_junction_idx = pipe["t_junction"]

    yp = getvariable(wm.model, :yp)[pipe_idx]
    yn = getvariable(wm.model, :yn)[pipe_idx]
    f = getvariable(wm.model, :f)[pipe_idx]

    max_flow = wm.data["max_flow"]
    hd_max = pipe["hd_max"]
    hd_min = pipe["hd_min"]
    w = pipe["resistance"]

    c1 = @constraint(wm.model, -(1-yp)*max_flow) <= f)
    c2 = @constraint(wm.model, f <= (1-yn)*max_flow)

    # c1 = @constraint(wm.model, -(1-yp)*min(max_flow, sqrt(w*max(hd_max, abs(hd_min)))) <= f)
    # c2 = @constraint(wm.model, f <= (1-yn)*min(max_flow, sqrt(w*max(hd_max, abs(hd_min)))))

    return Set([c1, c2])    # checked (though one is reveresed of the other)
end

#Bounds on pipe flow (by diameter)
function constraint_on_off_pipe_flow_bounds_diameter{T}(wm::GenericWaterModel{T}, pipe)
    pipe_idx = pipe["index"]
    i_junction_idx = pipe["f_junction"]
    j_junction_idx = pipe["t_junction"]
    D = pipe["diameter"]

    yp = getvariable(wm.model, :yp)[pipe_idx]
    yn = getvariable(wm.model, :yn)[pipe_idx]
    f = getvariable(wm.model, :f)[pipe_idx]

    max_flow = (pi/4)*(D^2)
    #max_flow = wm.data["max_flow"]
    hd_max = pipe["hd_max"]
    hd_min = pipe["hd_min"]
    # w = pipe["resistance"]

    c1 = @constraint(wm.model, -(1-yp)*max_flow) <= f)
    c2 = @constraint(wm.model, f <= (1-yn)*max_flow)

    # c1 = @constraint(wm.model, -(1-yp)*min(max_flow, sqrt(w*max(hd_max, abs(hd_min)))) <= f)
    # c2 = @constraint(wm.model, f <= (1-yn)*min(max_flow, sqrt(w*max(hd_max, abs(hd_min)))))

    return Set([c1, c2])    # checked (though one is reveresed of the other)
end


# standard flow balance equation where demand and production is fixed
function constraint_junction_flow_balance{T}(wm::GenericWaterModel{T}, junction)
    i = junction["index"]
    junction_branches = wm.set.junction_connections[i]

    f_branches = collect(keys(filter( (a, connection) -> connection["f_junction"] == i, wm.set.connections)))
    t_branches = collect(keys(filter( (a, connection) -> connection["t_junction"] == i, wm.set.connections)))

    h = getvariable(wm.model, :h_water)
    f = getvariable(wm.model, :f)

    c = @constraint(wm.model, junction["demand"] == sum(f[a] for a in t_branches) - sum(f[a] for a in f_branches) )

    return Set([c])
end

#Total head is lower bounded by elevation
function constraint_head_bound_elevation{T}(wm::GenericWaterModel{T}, junction)
    i = junction["index"]
    Hi = junction["elevation"]

    hi = getvariable(wm.model, :h_water)[i]

    c = @constraint(wm.model,hi >= Hi)

    return Set([c])
end


# # constraints on flow across valves
# function constraint_valve_flow{T}(wm::GenericWaterModel{T}, valve)
#     valve_idx = valve["index"]
#
#     f = getvariable(wm.model, :f)[valve_idx]
#     v = getvariable(wm.model, :valve)[valve_idx]
#
#     max_flow = wm.data["max_flow"]
#
#     c = @constraint(wm.model, -max_flow*v <= f <= max_flow*v)
#
#     return Set([c])
# end
#
# function constraint_valve_pressure_drop{T}(wm::GenericWaterModel{T})
#   hi = getvariable(wm.model, :h_water)[i_junction_idx]
#   hj = getvariable(wm.model, :h_water)[j_junction_idx]
#   v = getvariable(wm.model, :valve)[valve_idx]
#
#   hd_min = pipe["hd_min"]
#   hd_max = pipe["hd_max"]
#
#   c = @constraint(wm.model,(1-v)*hd_min <= hi - hj <=(1-v)*hd_max)
#
#   return Set([c])
# end
