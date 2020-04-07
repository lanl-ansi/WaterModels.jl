function run_owf(network, model_constructor, optimizer; kwargs...)
    return run_model(network, model_constructor, optimizer, post_owf; kwargs...)
end

function post_owf(wm::AbstractWaterModel)
    # Create head loss functions, if necessary.
    function_head_loss(wm)

    # Physical variables.
    variable_head(wm)
    variable_head_gain(wm)
    variable_flow(wm)
    variable_volume(wm)

    # Component-specific variables.
    variable_check_valve(wm)
    variable_pump_operation(wm)
    variable_reservoir(wm)
    variable_tank(wm)

    for (a, pipe) in ref(wm, :pipe)
        # TODO: Call this something other than status.
        if pipe["status"] == "CV"
            constraint_check_valve(wm, a)
            constraint_head_loss_check_valve(wm, a)
        else
            constraint_head_loss_pipe(wm, a)
        end
    end

    for a in ids(wm, :pump)
        constraint_head_gain_pump(wm, a)
    end

    # Flow conservation at all nodes.
    for (i, node) in ref(wm, :node)
        constraint_flow_conservation(wm, i)
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

    objective_owf(wm)
end

function run_mn_owf(file, model_constructor, optimizer; kwargs...)
    return run_model(file, model_constructor, optimizer, post_mn_owf; multinetwork=true, kwargs...)
end

function post_mn_owf(wm::AbstractWaterModel)
    # Create head loss functions, if necessary.
    function_head_loss(wm)

    for (n, network) in nws(wm)
        # Physical variables.
        variable_head(wm, nw=n)
        variable_head_gain(wm, nw=n)
        variable_flow(wm, nw=n)
        variable_volume(wm, nw=n)

        # Component-specific variables.
        variable_check_valve(wm, nw=n)
        variable_pump_operation(wm, nw=n)
        variable_reservoir(wm, nw=n)
        variable_tank(wm, nw=n)

        for (a, pipe) in ref(wm, :pipe, nw=n)
            # TODO: Call this something other than status.
            if pipe["status"] == "CV"
                constraint_check_valve(wm, a, nw=n)
                constraint_head_loss_check_valve(wm, a, nw=n)
            else
                constraint_head_loss_pipe(wm, a, nw=n)
            end
        end

        for a in ids(wm, :pump, nw=n)
            constraint_head_gain_pump(wm, a, nw=n)
        end

        # Flow conservation at all nodes.
        for (i, node) in ref(wm, :node, nw=n)
            constraint_flow_conservation(wm, i, nw=n)
        end

        # Set source node hydraulic heads.
        for (i, reservoir) in ref(wm, :reservoir, nw=n)
            constraint_source_head(wm, i, nw=n)
            constraint_source_flow(wm, i, nw=n)
        end

        for (i, junction) in ref(wm, :junction, nw=n)
            if junction["demand"] > 0.0
                constraint_sink_flow(wm, i, nw=n)
            elseif junction["demand"] < 0.0
                constraint_source_flow(wm, i, nw=n)
            end
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

    for n_2 in network_ids[2:end]
        # Set tank states after the initial time step.
        for i in ids(wm, :tank, nw=n_2)
            constraint_tank_state(wm, i, n_1, n_2)
        end

        n_1 = n_2
    end

    # Ensure tanks recover their initial volume.
    n_1, n_f = [network_ids[1], network_ids[end]]

    for i in ids(wm, :tank, nw=n_f)
        constraint_recover_volume(wm, i, n_1, n_f)
    end

    objective_owf(wm)
end
