function run_des(network, model_constructor, optimizer; kwargs...)
    return run_model(network, model_constructor, optimizer, build_des; kwargs...)
end

function build_des(wm::AbstractWaterModel)
    # Create head loss functions, if necessary.
    function_head_loss(wm)

    # Physical variables.
    variable_head(wm)
    variable_head_gain(wm)
    variable_flow(wm)

    # Component-specific variables.
    variable_reservoir(wm)

    # Add the network design objective.
    objective_des(wm)

    # Flow conservation at all nodes.
    for (i, node) in ref(wm, :node)
        constraint_flow_conservation(wm, i)
        constraint_node_directionality(wm, i)
    end

    # Head loss along fixed (non-design) pipes.
    for (a, pipe) in ref(wm, :pipe)
        constraint_pipe_head_loss(wm, a)
    end

    # Head loss along design pipes.
    for (a, pipe) in ref(wm, :des_pipe)
        constraint_pipe_head_loss_des(wm, a)
    end

    # Head loss along pipes with check valves.
    for (a, check_valve) in ref(wm, :check_valve)
        constraint_check_valve_head_loss(wm, a)
    end

    # Head loss along pipes with shutoff valves.
    for (a, shutoff_valve) in ref(wm, :shutoff_valve)
        constraint_shutoff_valve_head_loss(wm, a)
    end

    # Set source node hydraulic heads.
    for (i, reservoir) in ref(wm, :reservoir)
        constraint_source_directionality(wm, reservoir["node"])
    end

    # Constrain flow directions based on demand.
    for (i, junction) in ref(wm, :junction)
        if !junction["dispatchable"] && junction["demand"] > 0.0
            constraint_sink_directionality(wm, junction["node"])
        elseif junction["dispatchable"] && junction["demand_min"] > 0.0
            constraint_sink_directionality(wm, junction["node"])
        elseif !junction["dispatchable"] && junction["demand"] < 0.0
            constraint_source_directionality(wm, junction["node"])
        elseif junction["dispatchable"] && junction["demand_max"] < 0.0
            constraint_source_directionality(wm, junction["node"])
        end
    end
end
