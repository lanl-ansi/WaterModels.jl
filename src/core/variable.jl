########################################################################
# This file defines commonly-used variables for water systems models.
########################################################################

#function variable_head_ne(wm::GenericWaterModel, n::Int = wm.cnw)
#
#end

function variable_flow(wm::GenericWaterModel, n::Int = wm.cnw)
    lbs, ubs = calc_flow_bounds(wm.ref[:nw][n][:pipes])
    wm.var[:nw][n][:q] = @variable(wm.model, [id in keys(wm.ref[:nw][n][:pipes])],
                                   lowerbound = lbs[id], upperbound = ubs[id],
                                   basename = "q_$(n)", start = 1.0e-4)
end

function variable_head_common(wm::GenericWaterModel, n::Int = wm.cnw)
    lbs, ubs = calc_head_bounds(wm.ref[:nw][n][:junctions], wm.ref[:nw][n][:reservoirs])
    junction_ids = [key for key in keys(wm.ref[:nw][n][:junctions])]
    reservoir_ids = [key for key in keys(wm.ref[:nw][n][:reservoirs])]
    ids = [junction_ids; reservoir_ids]

    # Add the head variables to the model.
    wm.var[:nw][n][:h] = @variable(wm.model, [i in ids], lowerbound = lbs[i],
                                   upperbound = ubs[i], basename = "h_$(n)",
                                   start = ubs[i])
end

"Variables associated with building pipes."
function variable_pipe_ne_common(wm::GenericWaterModel, n::Int = wm.cnw)
    # Set up required data to initialize junction variables.
    pipe_ids = [key for key in keys(wm.ref[:nw][n][:ne_pipe])]
    wm.var[:nw][n][:psi] = Dict{Int, Any}()
    wm.var[:nw][n][:gamma] = Dict{Int, Any}()
    wm.var[:nw][n][:gamma_sum] = Dict{Int, Any}()

    for (pipe_id, pipe) in wm.ref[:nw][n][:ne_pipe]
        diameters = [d["diameter"] for d in pipe["diameters"]]

        # Create binary variables associated with whether or not a diameter is used.
        wm.var[:nw][n][:psi][pipe_id] = @variable(wm.model, [d in diameters], category = :Bin,
                                                  basename = "psi_$(n)_$(pipe_id)", start = 0)
        setvalue(wm.var[:nw][n][:psi][pipe_id][diameters[end]], 1)

        # Create a variable that corresponds to the selection of lambda.
        h_i = wm.var[:nw][n][:h][parse(Int, pipe["node1"])]
        h_j = wm.var[:nw][n][:h][parse(Int, pipe["node2"])]

        # Get additional data related to the variables.
        hij_lb = getlowerbound(h_i) - getupperbound(h_j)
        hij_ub = getupperbound(h_i) - getlowerbound(h_j)

        if pipe["flow_direction"] == POSITIVE
            hij_lb = max(0.0, hij_lb)
        elseif pipe["flow_direction"] == NEGATIVE
            hij_ub = min(0.0, hij_ub)
        end

        lbs = Dict(d => hij_lb / calc_friction_factor_hw_ne(pipe, d) for d in diameters)
        ubs = Dict(d => hij_ub / calc_friction_factor_hw_ne(pipe, d) for d in diameters)

        min_gamma, max_gamma = [minimum(collect(values(lbs))), maximum(collect(values(ubs)))]
        wm.var[:nw][n][:gamma][pipe_id] = @variable(wm.model, [d in diameters],
                                                    lowerbound = lbs[d], upperbound = ubs[d],
                                                    start = 0.5 * (ubs[d] + lbs[d]))
        wm.var[:nw][n][:gamma_sum][pipe_id] = @variable(wm.model, lowerbound = min_gamma,
                                                        upperbound = max_gamma,
                                                        start = 0.5 * (min_gamma + max_gamma))
    end
end

function variable_objective_ne(wm::GenericWaterModel, n::Int = wm.cnw)
    wm.var[:nw][n][:objective] = @variable(wm.model, basename = "objective_$(n)",
                                           lowerbound = 0.0, start = 0.0)
end
