function run_wf(network, model_constructor, optimizer; kwargs...)
    return run_model(network, model_constructor, optimizer, post_wf; kwargs...)
end

function post_wf(wm::AbstractWaterModel)
    # Create head loss functions, if necessary.
    function_head_loss(wm)

    # Physical variables.
    variable_head(wm, bounded=false)
    variable_flow(wm, bounded=false)
    variable_volume(wm, bounded=false)

    # Component-specific variables.
    variable_reservoir(wm)
    variable_tank(wm)

    # Head loss at pipes.
    for (a, pipe) in ref(wm, :pipe)
        constraint_head_loss_pipe(wm, a)
    end

    # Head gain along pumps (forced to be on).
    for (a, pump) in ref(wm, :pump)
        constraint_head_gain_pump(wm, a, force_on=true)
    end

    # Flow conservation at all nodes.
    for (i, node) in ref(wm, :node)
        constraint_flow_conservation(wm, i)
    end

    # Set source node hydraulic heads.
    for (i, reservoir) in ref(wm, :reservoir)
        constraint_source_head(wm, i)
        constraint_source_flow(wm, i)
    end

    for (i, junction) in ref(wm, :junction)
        # TODO: The conditional may be redundant, here.
        if junction["demand"] > 0.0
            constraint_sink_flow(wm, i)
        end
    end

    for i in ids(wm, :tank)
        # Link tank volume variables with tank head variables.
        constraint_link_volume(wm, i)

        # Set the initial tank volume.
        constraint_tank_state(wm, i)
    end

    # Add the objective.
    objective_wf(wm)
end

function run_mn_wf(file, model_constructor, optimizer; kwargs...)
    return run_model(file, model_constructor, optimizer, post_mn_wf; multinetwork=true, kwargs...)
end

function post_mn_wf(wm::AbstractWaterModel)
    # Create head loss functions, if necessary.
    function_head_loss(wm)

    for (n, network) in nws(wm)
        # Physical variables.
        variable_head(wm, nw=n, bounded=false)
        variable_flow(wm, nw=n, bounded=false)
        variable_volume(wm, nw=n, bounded=false)

        # Component-specific variables.
        variable_reservoir(wm, nw=n)
        variable_tank(wm, nw=n)

        # Head loss along pipes.
        for (a, pipe) in ref(wm, :pipe, nw=n)
            constraint_head_loss_pipe(wm, a, nw=n)
        end

        # Head gain along pumps (forced to be on).
        for (a, pump) in ref(wm, :pump, nw=n)
            constraint_head_gain_pump(wm, a, force_on=true, nw=n)
        end

        # Flow conservation at all nodes.
        for (i, node) in ref(wm, :node, nw=n)
            constraint_flow_conservation(wm, i, nw=n)
        end

        # Set source node hydraulic heads.
        for (i, reservoir) in ref(wm, :reservoir, nw=n)
            constraint_source_head(wm, i, nw=n)
        end

        # Link tank volume variables with tank head variables.
        for i in ids(wm, :tank, nw=n)
            constraint_link_volume(wm, i, nw=n)
        end
    end

    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(nw_ids(wm)))

    # Start with the first network, representing the initial time step.
    n_1 = network_ids[1]

    # Set initial conditions of tanks.
    for i in ids(wm, :tank, nw=n_1)
        constraint_tank_state(wm, i, nw=n_1)
    end

    # Set tank states after the initial time step.
    for n_2 in network_ids[2:end]
        for i in ids(wm, :tank, nw=n_2)
            constraint_tank_state(wm, i, n_1, n_2)
        end

        n_1 = n_2
    end

    # Add the objective.
    objective_wf(wm)
end
