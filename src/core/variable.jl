########################################################################
# This file defines commonly-used variables for water systems models.
########################################################################

function max_flow(ref, id)
    # Diameter assumes original units of millimeters.
    diameter = ref[:pipes][id]["diameter"] / 1000.0

    # A literature-based guess at v_max (meters per second).
    v_max = 10.0

    # Return the maximum flow.
    return 1.0e6 #(pi / 4.0) * v_max * diameter^2
end

function variable_flow{T}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    wm.var[:nw][n][:q] = @variable(wm.model, [id in keys(wm.ref[:nw][n][:pipes])],
                                   lowerbound = -max_flow(wm.ref[:nw][n], id),
                                   upperbound = max_flow(wm.ref[:nw][n], id),
                                   basename = "q_$(n)")

    wm.var[:nw][n][:gamma] = @variable(wm.model, [id in keys(wm.ref[:nw][n][:pipes])],
                                       lowerbound = -1000.0, upperbound = 1000.0,
                                       basename = "gamma_$(n)")
end

function variable_flow_direction{T}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    # Create variables that correspond to flow moving from i to j.
    wm.var[:nw][n][:yp] = @variable(wm.model, [id in keys(wm.ref[:nw][n][:pipes])],
                                    category = :Int, basename = "yp_$(n)",
                                    lowerbound = 0, upperbound = 1)

    # Create variables that correspond to flow moving from j to i.
    wm.var[:nw][n][:yn] = @variable(wm.model, [id in keys(wm.ref[:nw][n][:pipes])],
                                    category = :Int, basename = "yn_$(n)",
                                    lowerbound = 0, upperbound = 1)
end

function variable_head{T}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    # Set up required data to initialize junction variables.
    junction_ids = [key for key in keys(wm.ref[:nw][n][:junctions])]
    junction_lbs = Dict(id => wm.ref[:nw][n][:junctions][id]["elev"] for id in junction_ids)
    junction_ubs = Dict(id => 1000.0 for id in junction_ids)

    # Set up required data to initialize reservoir variables.
    reservoir_ids = [key for key in keys(wm.ref[:nw][n][:reservoirs])]
    reservoir_lbs = Dict(id => wm.ref[:nw][n][:reservoirs][id]["head"] for id in reservoir_ids)
    reservoir_ubs = Dict(id => wm.ref[:nw][n][:reservoirs][id]["head"] for id in reservoir_ids)

    # Create arrays comprising both types of components.
    ids = [junction_ids; reservoir_ids]
    lbs = merge(junction_lbs, reservoir_lbs)
    ubs = merge(junction_ubs, reservoir_ubs)

    # Add the head variables to the model.
    wm.var[:nw][n][:h] = @variable(wm.model, [i in ids], lowerbound = lbs[i],
                                   upperbound = ubs[i], basename = "h_$(n)")
end

#function variable_head{T}(wm::GenericWaterModel{T})
#    wm.var[:nw][n][:q] = @variable(wm.model, [i in keys(wm.ref[:junction])], basename="h", lowerbound=wm.ref[:junction][i]["head"], start = getstart(wm.ref[:junction], i, "h_start", wm.ref[:junction][i]["head"]))
#    # wm.var[:h] = @variable(wm.model, [i in keys(wm.ref[:junction])], basename="h", lowerbound=wm.ref[:junction][i]["elevation"]^2, upperbound=wm.ref[:junction][i]["pmax"]^2, start = getstart(wm.ref[:junction], i, "h_start", wm.ref[:junction][i]["pmin"]^2))
#    wm.var[:h] = @variable(wm.model, [i in keys(wm.ref[:junction])], basename="h", lowerbound=wm.ref[:junction][i]["head"], start = getstart(wm.ref[:junction], i, "h_start", wm.ref[:junction][i]["head"]))
#end

#"extracts the start value"
#function getstart(set, item_key, value_key, default = 0.0)
#    return get(get(set, item_key, Dict()), value_key, default)
#end
#
#" variables associated with head "
#function variable_head{T}(wm::GenericWaterModel{T})
#    # wm.var[:h] = @variable(wm.model, [i in keys(wm.ref[:junction])], basename="h", lowerbound=wm.ref[:junction][i]["elevation"]^2, upperbound=wm.ref[:junction][i]["pmax"]^2, start = getstart(wm.ref[:junction], i, "h_start", wm.ref[:junction][i]["pmin"]^2))
#    wm.var[:h] = @variable(wm.model, [i in keys(wm.ref[:junction])], basename="h", lowerbound=wm.ref[:junction][i]["head"], start = getstart(wm.ref[:junction], i, "h_start", wm.ref[:junction][i]["head"]))
#end
#
#" variables associated with flux "
#function variable_flux{T}(wm::GenericWaterModel{T})
#    max_flow = wm.ref[:max_flow]
#    wm.var[:f] = @variable(wm.model, [i in keys(wm.ref[:connection])], basename="f", lowerbound=-max_flow, upperbound=max_flow, start = getstart(wm.ref[:connection], i, "f_start", 0))
#end
#
## " variables associated with flux in expansion planning "
## function variable_flux_ne{T}(gm::GenericWaterModel{T})
##     max_flow = gm.ref[:max_flow]
##     gm.var[:f_ne] = @variable(gm.model, [i in keys(gm.ref[:ne_connection])], basename="f_ne", lowerbound=-max_flow, upperbound=max_flow, start = getstart(gm.ref[:ne_connection], i, "f_start", 0))
## end
#
#" variables associated with direction of flow on the connections "
#function variable_connection_direction{T}(wm::GenericWaterModel{T})
#    wm.var[:yp] = @variable(wm.model, [l in keys(wm.ref[:connection])], category = :Int, basename="yp", lowerbound=0, upperbound=1, start = getstart(wm.ref[:connection], l, "yp_start", 1.0))
#    wm.var[:yn] = @variable(wm.model, [l in keys(wm.ref[:connection])], category = :Int, basename="yn", lowerbound=0, upperbound=1, start = getstart(wm.ref[:connection], l, "yn_start", 0.0))
#end
#
## " variables associated with direction of flow on the connections "
## function variable_connection_direction_ne{T}(gm::GenericGasModel{T})
##      gm.var[:yp_ne] = @variable(gm.model, [l in keys(gm.ref[:ne_connection])], category = :Int,basename="yp_ne", lowerbound=0, upperbound=1, start = getstart(gm.ref[:ne_connection], l, "yp_start", 1.0))
##      gm.var[:yn_ne] = @variable(gm.model, [l in keys(gm.ref[:ne_connection])], category = :Int, basename="yn_ne", lowerbound=0, upperbound=1, start = getstart(gm.ref[:ne_connection], l, "yn_start", 0.0))
## end
#
## " variables associated with building pipes "
## function variable_pipe_ne{T}(gm::GenericGasModel{T})
##     gm.var[:zp] = @variable(gm.model, [l in keys(gm.ref[:ne_pipe])], category = :Int, basename="zp", lowerbound=0, upperbound=1, start = getstart(gm.ref[:ne_connection], l, "zp_start", 0.0))
## end
#
## " variables associated with building compressors "
## function variable_compressor_ne{T}(gm::GenericGasModel{T})
##     gm.var[:zc] = @variable(gm.model, [l in keys(gm.ref[:ne_compressor])], category = :Int,basename="zc", lowerbound=0, upperbound=1, start = getstart(gm.ref[:ne_connection], l, "zc_start", 0.0))
## end
#
#" 0-1 variables associated with operating valves "
#function variable_valve_operation{T}(wm::GenericWaterModel{T})
#    wm.var[:v] = @variable(wm.model, [l in [collect(keys(wm.ref[:valve])); collect(keys(wm.ref[:control_valve]))]], category = :Int, basename="v", lowerbound=0, upperbound=1, start = getstart(wm.ref[:connection], l, "v_start", 1.0))
#end
#
## " variables associated with demand "
## function variable_load{T}(gm::GenericGasModel{T})
##     load_set = filter(i -> gm.ref[:junction][i]["qlmin"] != gm.ref[:junction][i]["qlmax"], collect(keys(gm.ref[:junction])))
##     gm.var[:ql] = @variable(gm.model, [i in load_set], basename="ql", lowerbound=gm.ref[:junction][i]["qlmin"], upperbound=gm.ref[:junction][i]["qlmax"], start = getstart(gm.ref[:junction], i, "ql_start", 0.0))
## end
#
## " variables associated with production "
## function variable_production{T}(gm::GenericGasModel{T})
##     prod_set = filter(i -> gm.ref[:junction][i]["qgmin"] != gm.ref[:junction][i]["qgmax"], collect(keys(gm.ref[:junction])))
##     gm.var[:qg] = @variable(gm.model, [i in prod_set], basename="qg", lowerbound=gm.ref[:junction][i]["qgmin"], upperbound=gm.ref[:junction][i]["qgmax"], start = getstart(gm.ref[:junction], i, "qg_start", 0.0))
## end
