########################################################################
# This file defines commonly-used variables for water systems models.
########################################################################

function variable_head(wm::GenericWaterModel, n_n::Int=wm.cnw)
    # Get indices for all network nodes.
    nodes = [collect(ids(wm, n_n, :junctions));
             collect(ids(wm, n_n, :reservoirs))]

    # Get the bounds associated with heads at nodes.
    lbs, ubs = calc_head_bounds(wm, n_n)

    # Initialize variables associated with head.
    wm.var[:nw][n_n][:h] = @variable(wm.model, [i in nodes], lowerbound = lbs[i],
                                     upperbound = ubs[i], start = ubs[i],
                                     category = :Cont, basename = "h_$(n_n)")
end

function variable_head_difference(wm::GenericWaterModel, n_n::Int=wm.cnw)
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n_n, :connection))

    # Get the bounds for the head difference variables.
    lbs, ubs = calc_head_difference_bounds(wm, n_n)

    # Initialize variables associated with negative head differences.
    wm.var[:nw][n_n][:dhn] = @variable(wm.model, [a in arcs], lowerbound = 0.0,
                                       upperbound = abs(lbs[a]), start = abs(lbs[a]),
                                       category = :Cont, basename = "dhn_$(n_n)")

    # Initialize variables associated with positive head differences.
    wm.var[:nw][n_n][:dhp] = @variable(wm.model, [a in arcs], lowerbound = 0.0,
                                       upperbound = abs(ubs[a]), start = abs(ubs[a]),
                                       category = :Cont, basename = "dhp_$(n_n)")
end

function variable_flow_direction(wm::GenericWaterModel, n_n::Int=wm.cnw)
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n_n, :connection))

    # Initialize variables associated with flow direction. If this variable is
    # equal to one, the flow direction is from i to j. If it is equal to zero,
    # the flow direction is from j to i.
    wm.var[:nw][n_n][:dir] = @variable(wm.model, [a in arcs], start = 1,
                                       category = :Bin, basename = "dir_$(n_n)")

    # Fix direction bounds if they can be fixed.
    dh_lbs, dh_ubs = calc_head_difference_bounds(wm, n_n)
    for (a, connection) in wm.ref[:nw][n_n][:connection]
        if (dh_lbs[a] >= 0.0)
            setlowerbound(wm.var[:nw][n_n][:dir][a], 1)
            setvalue(wm.var[:nw][n_n][:dir][a], 1)
        elseif (dh_ubs[a] <= 0.0)
            setupperbound(wm.var[:nw][n_n][:dir][a], 0)
            setvalue(wm.var[:nw][n_n][:dir][a], 0)
        end
    end
end

function variable_resistance(wm::GenericWaterModel, n_n::Int=wm.cnw)
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n_n, :connection))

    # Initialize variables associated with flow direction. If this variable is
    # equal to one, the flow direction is from i to j. If it is equal to zero,
    # the flow direction is from j to i.
    wm.var[:nw][n_n][:xr] = Dict{Int, Array{Variable, 1}}()

    for (a, connection) in wm.ref[:nw][n_n][:connection]
        R_a = wm.ref[:nw][n_n][:resistance][a]

        wm.var[:nw][n_n][:xr][a] = @variable(wm.model, [r in 1:length(R_a)],
                                             start = 0, category = :Bin,
                                             basename = "xr_$(n_n)_$(a)")

        setvalue(wm.var[:nw][n_n][:xr][a][1], 1)
    end
end
