function solve_ne(network, model_constructor, optimizer; kwargs...)
    return solve_model(network, model_constructor, optimizer, build_ne; kwargs...)
end


function solve_mn_ne(file, model_constructor, optimizer; kwargs...)
    return solve_model(file, model_constructor, optimizer, build_mn_ne; multinetwork=true, kwargs...)
end


function build_ne(wm::AbstractWaterModel)
    # Build the water flow problem.
    build_wf(wm)

    # Add the optimal network expansion objective.
    objective_ne(wm)
end


function build_mn_ne(wm::AbstractWaterModel)
    # Build the water flow problem.
    build_mn_wf(wm)

    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(nw_ids(wm)))

    # Ensure tanks recover their initial volume.
    n_1, n_f = network_ids[1], network_ids[end]

    if(haskey(wm.ref[:it][wm_it_sym],:tank_volume_recovery_time_points))
        tank_volume_recovery_time_points = wm.ref[:it][wm_it_sym][:tank_volume_recovery_time_points]
    else
        tank_volume_recovery_time_points  = Set([])
    end
    union(tank_volume_recovery_time_points ,n_f)
    for n_tank in tank_volume_recovery_time_points
        for i in ids(wm, n_tank, :tank)
            constraint_tank_volume_recovery(wm, i, n_1, n_tank)
        end
    end

    # Add the optimal network expansion objective.
    objective_ne(wm)
end
