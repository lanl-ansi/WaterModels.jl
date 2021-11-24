function solve_owf(network, model_constructor, optimizer; kwargs...)
    return solve_model(network, model_constructor, optimizer, build_owf; kwargs...)
end


function solve_mn_owf(file, model_constructor, optimizer; kwargs...)
    return solve_model(file, model_constructor, optimizer, build_mn_owf; multinetwork=true, kwargs...)
end


function solve_mn_owf_switching(file, model_constructor, optimizer; kwargs...)
    return solve_model(file, model_constructor, optimizer, build_mn_owf_switching; multinetwork=true, kwargs...)
end


function build_owf(wm::AbstractWaterModel)
    # Build the water flow problem.
    build_wf(wm)

    # Add the optimal water flow objective.
    objective_owf(wm)
end


function build_mn_owf(wm::AbstractWaterModel)
    # Build the water flow problem.
    build_mn_wf(wm)

    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(nw_ids(wm)))

    # Ensure tanks recover their initial volume.
    n_1, n_f = network_ids[1], network_ids[end]

    for i in ids(wm, n_f, :tank)
        constraint_tank_volume_recovery(wm, i, n_1, n_f)
    end

    # Add the optimal water flow objective.
    objective_owf(wm)
end


function build_mn_owf_switching(wm::AbstractWaterModel)
    # Build the water flow problem.
    build_mn_wf_switching(wm)

    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(nw_ids(wm)))

    # Ensure tanks recover their initial volume.
    n_1, n_f = network_ids[1], network_ids[end]

    for i in ids(wm, n_f, :tank)
        constraint_tank_volume_recovery(wm, i, n_1, n_f)
    end

    # Add the optimal water flow objective.
    objective_owf(wm)
end
