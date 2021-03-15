function solve_owf(network, model_constructor, optimizer; kwargs...)
    return solve_model(network, model_constructor, optimizer, build_owf; kwargs...)
end


function run_owf(network, model_constructor, optimizer; kwargs...)
    Memento.warn(_LOGGER, "\"run_\" methods should be renamed \"solve_\" and will be deprecated in future versions.")
    return solve_owf(network, model_constructor, optimizer; kwargs...)
end


function solve_mn_owf(file, model_constructor, optimizer; kwargs...)
    return solve_model(file, model_constructor, optimizer, build_mn_owf; multinetwork=true, kwargs...)
end


function run_mn_owf(network, model_constructor, optimizer; kwargs...)
    Memento.warn(_LOGGER, "\"run_\" methods should be renamed \"solve_\" and will be deprecated in future versions.")
    return solve_mn_owf(network, model_constructor, optimizer; kwargs...)
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

function build_mn_owf_part(wm::AbstractWaterModel)
    # Build the water flow problem.
    build_mn_wf_decomp(wm)

    # Add the optimal water flow objective.
    #objective_owf(wm)
    objective_owf_decomp(wm)
end

function build_mn_owf_part_start(wm::AbstractWaterModel)
    # Build the water flow problem.
    build_mn_wf(wm)

    # Add the optimal water flow objective.
    objective_owf(wm)
end

function build_mn_owf_part_slack(wm::AbstractWaterModel)
    # Build the water flow problem.
    build_mn_wf_slack(wm)

    # Get all network IDs in the multinetwork.
    # network_ids = sort(collect(nw_ids(wm)))

    # n_1,n_f_2 = network_ids[1],network_ids[end]
    # n_2 = network_ids[1] + 1
    # n_f_1 = network_ids[end] - 1 

    # for i in ids(wm, n_f_2, :tank)
    #     constraint_tank_volume_with_slack_end(wm, i, n_f_1, n_f_2)
    # end

    # for i in ids(wm, n_2, :tank)
    #     constraint_tank_volume_with_slack_start(wm, i, n_1, n_2)
    # end

    # Add the optimal water flow objective.
    objective_owf(wm)
end

