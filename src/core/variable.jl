########################################################################
# This file defines commonly-used variables for water systems models.
########################################################################

function variable_flow{T}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    flow_min, flow_max = calc_flow_bounds(wm.ref[:nw][n][:pipes])
    wm.var[:nw][n][:q] = @variable(wm.model,
                                   [id in keys(wm.ref[:nw][n][:pipes])],
                                   lowerbound = flow_min[id],
                                   upperbound = flow_max[id],
                                   basename = "q_$(n)")

end

function variable_head_difference{T}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    diff_min, diff_max = calc_head_difference_bounds(wm.ref[:nw][n][:pipes])
    wm.var[:nw][n][:gamma] = @variable(wm.model,
                                       [id in keys(wm.ref[:nw][n][:pipes])],
                                       lowerbound = diff_min[id],
                                       upperbound = diff_max[id],
                                       basename = "gamma_$(n)")
end

function variable_flow_direction{T}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    # Create variables that correspond to flow moving from i to j.
    wm.var[:nw][n][:yp] = @variable(wm.model,
                                    [id in keys(wm.ref[:nw][n][:pipes])],
                                    lowerbound = 0,
                                    upperbound = 1,
                                    category = :Int,
                                    basename = "yp_$(n)")

    # Create variables that correspond to flow moving from j to i.
    wm.var[:nw][n][:yn] = @variable(wm.model,
                                    [id in keys(wm.ref[:nw][n][:pipes])],
                                    lowerbound = 0,
                                    upperbound = 1,
                                    category = :Int,
                                    basename = "yn_$(n)")
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
