##########################################################################################################
# The purpose of this file is to define commonly used and created variables used in water flow models
##########################################################################################################

# extracts the start value from data,
function getstart(set, item_key, value_key, default = 0.0)
    try
        return set[item_key][value_key]
    catch
        return default
    end
end

# variables associated with water head
function variable_water_head{T}(wm::GenericWaterModel{T})
    @variable(wm.model, wm.set.junctions[i]["hmin"]^2 <= h_water[i in gm.set.junction_indexes] <= wm.set.junctions[i]["hmax"]^2, start = getstart(wm.set.junctions, i, "p_start", wm.set.junctions[i]["pmin"]^2))
    return h_water
end

# variables associated with flux
function variable_flux{T}(wm::GenericWaterModel{T})
    max_flow = wm.data["max_flow"]
    @variable(wm.model, -max_flow <= f[i in wm.set.connection_indexes] <= max_flow, start = getstart(wm.set.connections, i, "f_start", 0))
    return f
end

# # variables associated with flux in expansion planning
# function variable_flux_ne{T}(wm::GenericWaterModel{T})
#     max_flow = wm.data["max_flow"]
#     @variable(wm.model, -max_flow <= f_ne[i in wm.set.new_connection_indexes] <= max_flow, start = getstart(wm.set.new_connections, i, "f_start", 0))
#     return f_ne
# end

# variables associated with direction of flow on the connections
function variable_connection_direction{T}(wm::GenericWaterModel{T})
    @variable(wm.model, 0 <= yp[l in wm.set.connection_indexes] <= 1, Int, start = getstart(wm.set.connections, l, "yp_start", 1.0))
    @variable(wm.model, 0 <= yn[l in wm.set.connection_indexes] <= 1, Int, start = getstart(wm.set.connections, l, "yn_start", 0.0))
    return yp, yn
end

# # variables associated with direction of flow on the connections
# function variable_connection_direction_ne{T}(wm::GenericWaterModel{T})
#     @variable(wm.model, 0 <= yp_ne[l in wm.set.new_connection_indexes] <= 1, Int, start = getstart(wm.set.new_connections, l, "yp_start", 1.0))
#     @variable(wm.model, 0 <= yn_ne[l in wm.set.new_connection_indexes] <= 1, Int, start = getstart(wm.set.new_connections, l, "yn_start", 0.0))
#     return yp_ne, yn_ne
# end

# # variables associated with building pipes
# function variable_pipe_ne{T}(wm::GenericWaterModel{T})
#     @variable(wm.model, 0 <= zp[l in wm.set.new_pipe_indexes] <= 1, Int, start = getstart(wm.set.new_connections, l, "zp_start", 0.0))
#     return zp
# end

# # variables associated with building compressors
# function variable_compressor_ne{T}(gm::GenericGasModel{T})
#     @variable(gm.model, 0 <= zc[l in gm.set.new_compressor_indexes] <= 1, Int, start = getstart(gm.set.new_connections, l, "zc_start", 0.0))
#     return zc
# end



# 0-1 variables associated with operating valves
function variable_valve_operation{T}(wm::GenericWaterModel{T})
    @variable(wm.model, 0 <= valve[l in wm.set.valve_indexes] <= 1, Int, start = getstart(wm.set.connections, l, "v_start", 1.0))
    return valve
end

# # variables associated with demand
# function variable_load{T}(gm::GenericGasModel{T})
#     load_set = filter(i -> gm.set.junctions[i]["qlmin"] != gm.set.junctions[i]["qlmax"], gm.set.junction_indexes)
#     @variable(gm.model, gm.set.junctions[i]["qlmin"] <= ql_gas[i in load_set] <= gm.set.junctions[i]["qlmax"], start = getstart(gm.set.junctions, i, "ql_start", 0.0))
#     return ql_gas
# end
#
# # variables associated with production
# function variable_production{T}(gm::GenericGasModel{T})
#     prod_set = filter(i -> gm.set.junctions[i]["qgmin"] != gm.set.junctions[i]["qgmax"], gm.set.junction_indexes)
#     @variable(gm.model, gm.set.junctions[i]["qgmin"] <= qg_gas[i in prod_set] <= gm.set.junctions[i]["qgmax"], start = getstart(gm.set.junctions, i, "qg_start", 0.0))
#     return qg_gas
# end
