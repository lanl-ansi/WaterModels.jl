function solve_des(network, model_constructor, optimizer; kwargs...)
    return solve_model(network, model_constructor, optimizer, build_des; kwargs...)
end

function build_des(wm::AbstractWaterModel)
    # Create head loss functions, if necessary.
    function_head_loss(wm)

    # Physical variables.
    variable_head(wm)
    variable_head_gain(wm)
    variable_flow(wm)
    variable_flow_des(wm)
    variable_volume(wm)

    # Component-specific variables.
    variable_check_valve(wm)
    variable_pump_operation(wm)
    variable_reservoir(wm)
    variable_tank(wm)

    # Add the network design objective.
    objective_des(wm)

    for (a, pipe) in ref(wm, :pipe_fixed)
        # TODO: Call this something other than status.
        if pipe["status"] == "CV"
            constraint_check_valve(wm, a)
            constraint_head_loss_check_valve(wm, a)
        else
            constraint_head_loss_pipe(wm, a)
        end
    end

    for (a, pipe) in ref(wm, :pipe_des)
        constraint_head_loss_pipe_des(wm, a)
    end

    for a in ids(wm, :pump)
        constraint_head_gain_pump(wm, a)
        constraint_pump_control(wm, a)
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

    # Add the energy conservation (primal-dual) constraint.
    constraint_energy_conservation(wm)
end
