########################################################################
# This file defines commonly-used variables for water systems models.
########################################################################

function variable_flow{T}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    lbs, ubs = calc_flow_bounds(wm.ref[:nw][n][:pipes])
    wm.var[:nw][n][:q] = @variable(wm.model, [id in keys(wm.ref[:nw][n][:pipes])],
                                   lowerbound = lbs[id], upperbound = ubs[id],
                                   basename = "q_$(n)", start = ubs[id])

end

function variable_head_common{T}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    # Set up required data to initialize junction variables.
    junction_ids = [key for key in keys(wm.ref[:nw][n][:junctions])]

    # Set up required data to initialize reservoir variables.
    reservoirs = wm.ref[:nw][n][:reservoirs]
    reservoir_ids = [key for key in keys(reservoirs)]
    reservoir_lbs = Dict(id => reservoirs[id]["head"] for id in reservoir_ids)
    reservoir_ubs = Dict(id => reservoirs[id]["head"] for id in reservoir_ids)

    # Set the elevation bounds (for junctions).
    # TODO: Increase the upper bound when pumps are in the system.
    junctions = wm.ref[:nw][n][:junctions]
    max_elev = maximum([junc["elev"] for junc in values(junctions)])
    max_head = maximum([res["head"] for res in values(reservoirs)])
    junction_lbs = Dict(junc["id"] => junc["elev"] for junc in values(junctions))
    junction_ubs = Dict(id => max(max_elev, max_head) for id in junction_ids)

    # Create arrays comprising both types of components.
    ids = [junction_ids; reservoir_ids]
    lbs = merge(junction_lbs, reservoir_lbs)
    ubs = merge(junction_ubs, reservoir_ubs)

    # Add the head variables to the model.
    wm.var[:nw][n][:h] = @variable(wm.model, [i in ids], lowerbound = lbs[i],
                                   upperbound = ubs[i], basename = "h_$(n)",
                                   start = ubs[i])
end

"Variables associated with building pipes."
function variable_pipe_ne{T}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    # Set up required data to initialize junction variables.
    pipe_ids = [key for key in keys(wm.ref[:nw][n][:ne_pipe])]
    wm.var[:nw][n][:psi] = Dict{String, Any}()
    wm.var[:nw][n][:gamma] = Dict{String, Any}()

    for (pipe_id, pipe) in wm.ref[:nw][n][:ne_pipe]
        diameters = [d["diameter"] for d in pipe["diameters"]]

        # Create binary variables associated with whether or not a diameter is used.
        wm.var[:nw][n][:psi][pipe_id] = @variable(wm.model, [d in diameters], category = :Bin,
                                                  basename = "psi_$(n)_$(pipe_id)", start = 0)
        setvalue(wm.var[:nw][n][:psi][pipe_id][diameters[end]], 1)

        # Create a variable that corresponds to the selection of lambda.
        lambdas = Dict(d => calc_friction_factor_hw_ne(pipe, d) for d in diameters)
        wm.var[:nw][n][:gamma][pipe_id] = @variable(wm.model, [d in diameters],
                                                    lowerbound = -1000.0 / lambdas[d],
                                                    upperbound = 1000.0 / lambdas[d],
                                                    start = 0.0)
    end
end
