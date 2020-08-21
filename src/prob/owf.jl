function run_owf(network, model_constructor, optimizer; kwargs...)
    return run_model(network, model_constructor, optimizer, build_owf; kwargs...)
end

function build_owf(wm::AbstractWaterModel)
    # Build the water flow problem.
    build_wf(wm)

    # Add the optimal water flow objective.
    objective_owf(wm)
end

function run_mn_owf(file, model_constructor, optimizer; kwargs...)
    return run_model(file, model_constructor, optimizer, build_mn_owf; multinetwork=true, kwargs...)
end

function build_mn_owf(wm::AbstractWaterModel)
    # Build the water flow problem.
    build_mn_wf(wm)

    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(nw_ids(wm)))

    # Ensure tanks recover their initial volume.
    n_1, n_f = network_ids[1], network_ids[end]

    for i in ids(wm, :tank; nw=n_f)
        constraint_recover_volume(wm, i, n_1, n_f)
    end

    # Add the optimal water flow objective.
    objective_owf(wm)
end
