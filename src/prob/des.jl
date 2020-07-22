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
    variable_flow_des(wm)

    # Component-specific variables.
    variable_reservoir(wm)

    # Add the network design objective.
    objective_des(wm)

    for (a, pipe) in ref(wm, :pipe_fixed)
        constraint_pipe_head_loss(wm, a)
    end

    for (a, check_valve) in ref(wm, :check_valve)
        constraint_check_valve_head_loss(wm, a)
    end

    for (a, shutoff_valve) in ref(wm, :shutoff_valve)
        constraint_shutoff_valve_head_loss(wm, a)
    end

    # Head loss along design pipes.
    for (a, pipe) in ref(wm, :pipe_des)
        constraint_pipe_head_loss_des(wm, a)
    end

    # Flow conservation at all nodes.
    for (i, node) in ref(wm, :node)
        constraint_flow_conservation(wm, i)
    end

    # Set source node hydraulic heads.
    for (i, reservoir) in ref(wm, :reservoir)
        constraint_reservoir_head(wm, i)
        constraint_source_flow(wm, i)
    end

    # Constrain flow directions based on demand.
    for (i, junction) in ref(wm, :junction)
        if junction["demand"] > 0.0
            constraint_sink_flow(wm, i)
        elseif junction["demand"] < 0.0
            constraint_source_flow(wm, i)
        end
    end
end
