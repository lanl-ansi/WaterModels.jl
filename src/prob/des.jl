function run_des(network, model_constructor, optimizer; kwargs...)
    return run_model(network, model_constructor, optimizer, build_des; kwargs...)
end


function build_des(wm::AbstractWaterModel)
    # Create head loss functions, if necessary.
    _function_head_loss(wm)

    # Physical variables.
    variable_head(wm)
    variable_flow(wm)
    variable_pump_head_gain(wm)

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
        constraint_pipe_flow(wm, a)
        constraint_pipe_head(wm, a)
        constraint_pipe_head_loss(wm, a)
    end

    # Constraints on design pipe flows, heads, and physics.
    for (a, des_pipe) in ref(wm, :des_pipe)
        constraint_on_off_des_pipe_flow(wm, a)
        constraint_on_off_des_pipe_head(wm, a)
        constraint_on_off_des_pipe_head_loss(wm, a)
    end

    # Selection of design pipes along unique arcs.
    for arc in ref(wm, :des_pipe_arcs)
        constraint_des_pipe_flow(wm, arc[1], arc[2])
        constraint_des_pipe_head(wm, arc[1], arc[2])
        constraint_des_pipe_selection(wm, arc[1], arc[2])
    end

    # Constraints on pump flows, heads, and physics.
    for (a, pump) in ref(wm, :pump)
        constraint_on_off_pump_head(wm, a)
        constraint_on_off_pump_head_gain(wm, a)
        constraint_on_off_pump_flow(wm, a)
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

    # Constraints on valve flows and heads.
    for (a, valve) in ref(wm, :valve)
        constraint_on_off_valve_head(wm, a)
        constraint_on_off_valve_flow(wm, a)
    end

    # Constraints on tank volumes.
    for i in ids(wm, :tank)
        # Set the initial tank volume.
        constraint_tank_volume(wm, i)
    end

    # Add the network design objective.
    objective_des(wm)
end
