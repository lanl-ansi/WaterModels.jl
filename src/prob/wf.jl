function solve_wf(network, model_constructor, optimizer; kwargs...)
    return solve_model(network, model_constructor, optimizer, build_wf; kwargs...)
end


function run_wf(network, model_constructor, optimizer; kwargs...)
    Memento.warn(_LOGGER, "\"run_\" methods should be renamed \"solve_\" and will be deprecated in future versions.")
    return solve_wf(network, model_constructor, optimizer; kwargs...)
end


function solve_mn_wf(file, model_constructor, optimizer; kwargs...)
    return solve_model(file, model_constructor, optimizer, build_mn_wf; multinetwork=true, kwargs...)
end


function run_mn_wf(network, model_constructor, optimizer; kwargs...)
    Memento.warn(_LOGGER, "\"run_\" methods should be renamed \"solve_\" and will be deprecated in future versions.")
    return solve_mn_wf(network, model_constructor, optimizer; kwargs...)
end


function build_wf(wm::AbstractWaterModel)
    # Create head loss functions, if necessary.
    _function_head_loss(wm)

    # Physical variables.
    variable_head(wm)
    variable_flow(wm)
    variable_pump_head_gain(wm)
    variable_pump_power(wm)

    # Indicator (status) variables.
    variable_des_pipe_indicator(wm)
    variable_pump_indicator(wm)
    variable_regulator_indicator(wm)
    variable_valve_indicator(wm)

    # Create flow-related variables for node attachments.
    variable_demand_flow(wm)
    variable_reservoir_flow(wm)
    variable_tank_flow(wm)

    # Flow conservation at all nodes.
    for (i, node) in ref(wm, :node)
        constraint_flow_conservation(wm, i)
        constraint_node_directionality(wm, i)
    end

    # Constraints on pipe flows, heads, and physics.
    for (a, pipe) in ref(wm, :pipe)
        constraint_pipe_head(wm, a)
        constraint_pipe_head_loss(wm, a)
        constraint_pipe_flow(wm, a)
    end

    # Selection of design pipes along unique arcs.
    for (k, arc) in ref(wm, :des_pipe_arc)
        constraint_des_pipe_flow(wm, k, arc[1], arc[2])
        constraint_des_pipe_head(wm, k, arc[1], arc[2])
        constraint_des_pipe_selection(wm, k, arc[1], arc[2])
    end

    # Constraints on design pipe flows, heads, and physics.
    for (a, des_pipe) in ref(wm, :des_pipe)
        constraint_on_off_des_pipe_head(wm, a)
        constraint_on_off_des_pipe_head_loss(wm, a)
        constraint_on_off_des_pipe_flow(wm, a)
    end

    # Constraints on pump flows, heads, and physics.
    for (a, pump) in ref(wm, :pump)
        constraint_on_off_pump_head(wm, a)
        constraint_on_off_pump_head_gain(wm, a)
        constraint_on_off_pump_flow(wm, a)
        constraint_on_off_pump_power(wm, a)
    end

    # Constraints on short pipe flows and heads.
    for (a, regulator) in ref(wm, :regulator)
        constraint_on_off_regulator_head(wm, a)
        constraint_on_off_regulator_flow(wm, a)
    end

    # Constraints on short pipe flows and heads.
    for (a, short_pipe) in ref(wm, :short_pipe)
        constraint_short_pipe_head(wm, a)
        constraint_short_pipe_flow(wm, a)
    end

    # Constraints on tank volumes.
    for (i, tank) in ref(wm, :tank)
        # Set the initial tank volume.
        constraint_tank_volume(wm, i)
    end

    # Constraints on valve flows and heads.
    for (a, valve) in ref(wm, :valve)
        constraint_on_off_valve_head(wm, a)
        constraint_on_off_valve_flow(wm, a)
    end

    # Add the objective.
    objective_wf(wm)
end


function build_mn_wf(wm::AbstractWaterModel)
    # Create head loss functions, if necessary.
    _function_head_loss(wm)

    for (n, network) in nws(wm)
        # Physical variables.
        variable_head(wm; nw=n)
        variable_flow(wm; nw=n)
        variable_pump_head_gain(wm; nw=n)
        variable_pump_power(wm; nw=n)

        # Indicator (status) variables.
        variable_des_pipe_indicator(wm; nw=n)
        variable_pump_indicator(wm; nw=n)
        variable_regulator_indicator(wm; nw=n)
        variable_valve_indicator(wm; nw=n)

        # Create flow-related variables for node attachments.
        variable_demand_flow(wm; nw=n)
        variable_reservoir_flow(wm; nw=n)
        variable_tank_flow(wm; nw=n)

        # Flow conservation at all nodes.
        for (i, node) in ref(wm, :node; nw=n)
            constraint_flow_conservation(wm, i; nw=n)
            constraint_node_directionality(wm, i; nw=n)
        end

        # Constraints on pipe flows, heads, and physics.
        for (a, pipe) in ref(wm, :pipe; nw=n)
            constraint_pipe_flow(wm, a; nw=n)
            constraint_pipe_head(wm, a; nw=n)
            constraint_pipe_head_loss(wm, a; nw=n)
        end

        # Constraints on design pipe flows, heads, and physics.
        for (a, des_pipe) in ref(wm, :des_pipe; nw=n)
            constraint_on_off_des_pipe_flow(wm, a; nw=n)
            constraint_on_off_des_pipe_head(wm, a; nw=n)
            constraint_on_off_des_pipe_head_loss(wm, a; nw=n)
        end

        # Selection of design pipes along unique arcs.
        for (k, arc) in ref(wm, :des_pipe_arc; nw=n)
            constraint_des_pipe_flow(wm, k, arc[1], arc[2]; nw=n)
            constraint_des_pipe_head(wm, k, arc[1], arc[2]; nw=n)
            constraint_des_pipe_selection(wm, k, arc[1], arc[2]; nw=n)
        end

        # Constraints on pump flows, heads, and physics.
        for (a, pump) in ref(wm, :pump; nw=n)
            constraint_on_off_pump_head(wm, a; nw=n)
            constraint_on_off_pump_head_gain(wm, a; nw=n)
            constraint_on_off_pump_flow(wm, a; nw=n)
            constraint_on_off_pump_power(wm, a; nw=n)
        end

        # Constraints on short pipe flows and heads.
        for (a, regulator) in ref(wm, :regulator; nw=n)
            constraint_on_off_regulator_head(wm, a; nw=n)
            constraint_on_off_regulator_flow(wm, a; nw=n)
        end

        # Constraints on short pipe flows and heads.
        for (a, short_pipe) in ref(wm, :short_pipe; nw=n)
            constraint_short_pipe_head(wm, a; nw=n)
            constraint_short_pipe_flow(wm, a; nw=n)
        end

        # Constraints on valve flows and heads.
        for (a, valve) in ref(wm, :valve; nw=n)
            constraint_on_off_valve_head(wm, a; nw=n)
            constraint_on_off_valve_flow(wm, a; nw=n)
        end
    end

    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(nw_ids(wm)))

    # Start with the first network, representing the initial time step.
    n_1 = network_ids[1]

    # Constraints on tank volumes.
    for (i, tank) in ref(wm, :tank; nw = n_1)
        # Set initial conditions of tanks.
        constraint_tank_volume(wm, i; nw = n_1)
    end

    # Constraints on tank volumes.
    for n_2 in network_ids[2:end]
        # Constrain tank volumes after the initial time step.
        for (i, tank) in ref(wm, :tank; nw = n_2)
            constraint_tank_volume(wm, i, n_1, n_2)
        end

        n_1 = n_2 # Update the first network used for integration.
    end

    # Add the objective.
    objective_wf(wm)
end
