function run_wf(network, model_constructor, optimizer; kwargs...)
    return run_model(network, model_constructor, optimizer, build_wf; kwargs...)
end

function build_wf(wm::AbstractWaterModel)
    # Create head loss functions, if necessary.
    function_head_loss(wm)

    # Physical variables.
    variable_head(wm)
    variable_head_gain(wm)
    variable_flow(wm)

    # Indicator (status) variables.
    variable_check_valve_indicator(wm)
    variable_pump_indicator(wm)
    variable_pressure_reducing_valve_indicator(wm)
    variable_shutoff_valve_indicator(wm)

    # Component-specific variables.
    variable_pump_operation(wm)
    variable_reservoir(wm)
    variable_tank(wm)

    # Flow conservation at all nodes.
    for (i, node) in ref(wm, :node)
        constraint_flow_conservation(wm, i)
    end

    # Head loss along pipes without valves.
    for (a, pipe) in ref(wm, :pipe)
        constraint_pipe_head_loss(wm, a)
    end

    # Head loss along pipes with check valves.
    for (a, check_valve) in ref(wm, :check_valve)
        constraint_check_valve_head_loss(wm, a)
    end

    # Head loss along pipes with shutoff valves.
    for (a, shutoff_valve) in ref(wm, :shutoff_valve)
        constraint_shutoff_valve_head_loss(wm, a)
    end

    # Head loss along pressure reducing valves.
    for a in ids(wm, :pressure_reducing_valve)
        constraint_prv_head_loss(wm, a)
    end

    # Head gain along pumps.
    for a in ids(wm, :pump)
        constraint_pump_head_gain(wm, a)
    end

    # Set source node hydraulic heads.
    for (i, reservoir) in ref(wm, :reservoir)
        constraint_reservoir_head(wm, i)
        constraint_source_directionality(wm, reservoir["node"])
    end

    # Constrain flow directions based on demand.
    for (i, junction) in ref(wm, :junction)
        if junction["demand"] > 0.0
            constraint_sink_directionality(wm, junction["node"])
        elseif junction["demand"] < 0.0
            constraint_source_directionality(wm, junction["node"])
        end
    end

    for i in ids(wm, :tank)
        # Set the initial tank volume.
        constraint_tank_state(wm, i)
    end

    # Add the objective.
    objective_wf(wm)
end

function run_mn_wf(file, model_constructor, optimizer; kwargs...)
    return run_model(file, model_constructor, optimizer, build_mn_wf; multinetwork=true, kwargs...)
end

function build_mn_wf(wm::AbstractWaterModel)
    # Create head loss functions, if necessary.
    function_head_loss(wm)

    for (n, network) in nws(wm)
        # Physical variables.
        variable_head(wm; nw=n)
        variable_head_gain(wm; nw=n)
        variable_flow(wm; nw=n)

        # Indicator (status) variables.
        variable_check_valve_indicator(wm; nw=n)
        variable_shutoff_valve_indicator(wm; nw=n)
        variable_pressure_reducing_valve_indicator(wm; nw=n)
        variable_pump_indicator(wm; nw=n)

        # Component-specific variables.
        variable_pump_operation(wm; nw=n)
        variable_reservoir(wm; nw=n)
        variable_tank(wm; nw=n)

        # Flow conservation at all nodes.
        for i in ids(wm, :node; nw=n)
            constraint_flow_conservation(wm, i; nw=n)
        end

        # Head loss along pipes without valves.
        for (a, pipe) in ref(wm, :pipe, nw=n)
            constraint_pipe_head_loss(wm, a, nw=n)
        end

        # Head loss along pipes with check valves.
        for (a, check_valve) in ref(wm, :check_valve, nw=n)
            constraint_check_valve_head_loss(wm, a, nw=n)
        end

        # Head loss along pipes with shutoff valves.
        for (a, shutoff_valve) in ref(wm, :shutoff_valve, nw=n)
            constraint_shutoff_valve_head_loss(wm, a, nw=n)
        end

        # Head loss along pressure reducing valves.
        for a in ids(wm, :pressure_reducing_valve; nw=n)
            constraint_prv_head_loss(wm, a; nw=n)
        end

        # Head gain along pumps.
        for a in ids(wm, :pump; nw=n)
            constraint_pump_head_gain(wm, a; nw=n)
        end

        # Constrain source node hydraulic heads and flow directions.
        for (i, reservoir) in ref(wm, :reservoir; nw=n)
            constraint_reservoir_head(wm, i; nw=n)
            constraint_source_directionality(wm, reservoir["node"]; nw=n)
        end

        # Constrain flow directions based on demand.
        for (i, junction) in ref(wm, :junction; nw=n)
            if junction["demand"] > 0.0
                constraint_sink_directionality(wm, junction["node"]; nw=n)
            elseif junction["demand"] < 0.0
                constraint_source_directionality(wm, junction["node"]; nw=n)
            end
        end
    end

    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(nw_ids(wm)))

    # Start with the first network, representing the initial time step.
    n_1 = network_ids[1]

    # Set initial conditions of tanks.
    for i in ids(wm, :tank; nw=n_1)
        constraint_tank_state(wm, i; nw=n_1)
    end

    for n_2 in network_ids[2:end]
        # Set tank states after the initial time step.
        for i in ids(wm, :tank; nw=n_2)
            constraint_tank_state(wm, i, n_1, n_2)
        end

        n_1 = n_2 # Update the first network used for integration.
    end

    # Add the objective.
    objective_wf(wm)
end
