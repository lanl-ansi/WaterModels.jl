function run_wf(network, model_constructor, optimizer; kwargs...)
    return run_model(network, model_constructor, optimizer, build_wf; kwargs...)
end

function build_wf(wm::AbstractWaterModel)
    # Create head loss functions, if necessary.
    function_head_loss(wm)

    # Physical variables.
    variable_head(wm)
    variable_flow(wm)
    variable_pump_head_gain(wm)

    # Indicator (status) variables.
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
        constraint_pipe_flow(wm, a)
        constraint_pipe_head_loss(wm, a)
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

    ## Head gain along pumps.
    #for a in ids(wm, :pump)
    #    constraint_pump_head_gain(wm, a)
    #end

    ## Set source node hydraulic heads.
    #for (i, reservoir) in ref(wm, :reservoir)
    #    constraint_reservoir_head(wm, i)
    #    constraint_source_directionality(wm, reservoir["node"])
    #end

    ## Constrain flow directions based on demand.
    #for (i, demand) in ref(wm, :demand)
    #    if !demand["dispatchable"] && demand["flow_rate"] > 0.0
    #        constraint_sink_directionality(wm, demand["node"])
    #    elseif demand["dispatchable"] && demand["demand_min"] > 0.0
    #        constraint_sink_directionality(wm, demand["node"])
    #    elseif !demand["dispatchable"] && demand["flow_rate"] < 0.0
    #        constraint_source_directionality(wm, demand["node"])
    #    elseif demand["dispatchable"] && demand["demand_max"] < 0.0
    #        constraint_source_directionality(wm, demand["node"])
    #    end
    #end

    #for i in ids(wm, :tank)
    #    # Set the initial tank volume.
    #    constraint_tank_state(wm, i)
    #end

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
        variable_flow(wm; nw=n)
        variable_pump_head_gain(wm; nw=n)

        # Indicator (status) variables.
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
            constraint_pipe_head(wm, a; nw=n)
            constraint_pipe_flow(wm, a; nw=n)
            constraint_pipe_head_loss(wm, a; nw=n)
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

        ## Head gain along pumps.
        #for a in ids(wm, :pump; nw=n)
        #    constraint_pump_head_gain(wm, a; nw=n)
        #end

        ## Constrain source node hydraulic heads and flow directions.
        #for (i, reservoir) in ref(wm, :reservoir; nw=n)
        #    constraint_reservoir_head(wm, i; nw=n)
        #    constraint_source_directionality(wm, reservoir["node"]; nw=n)
        #end

        ## Constrain flow directions based on demand.
        #for (i, demand) in ref(wm, :demand; nw=n)
        #    if !demand["dispatchable"] && demand["flow_rate"] > 0.0
        #        constraint_sink_directionality(wm, demand["node"]; nw=n)
        #    elseif !demand["dispatchable"] && demand["flow_rate"] < 0.0
        #        constraint_source_directionality(wm, demand["node"]; nw=n)
        #    end
        #end
    end

    ## Get all network IDs in the multinetwork.
    #network_ids = sort(collect(nw_ids(wm)))

    ## Start with the first network, representing the initial time step.
    #n_1 = network_ids[1]

    ## Set initial conditions of tanks.
    #for i in ids(wm, :tank; nw=n_1)
    #    constraint_tank_state(wm, i; nw=n_1)
    #end

    #for n_2 in network_ids[2:end]
    #    # Set tank states after the initial time step.
    #    for i in ids(wm, :tank; nw=n_2)
    #        constraint_tank_state(wm, i, n_1, n_2)
    #    end

    #    n_1 = n_2 # Update the first network used for integration.
    #end

    # Add the objective.
    objective_wf(wm)
end
