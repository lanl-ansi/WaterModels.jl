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

    for i in ids(wm, n_f, :tank)
        constraint_tank_volume_recovery(wm, i, n_1, n_f)
    end

    # Add the optimal network expansion objective.
    objective_ne(wm)
end