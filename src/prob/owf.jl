function run_owf(network, model_constructor, optimizer; kwargs...)
    return run_model(network, model_constructor, optimizer, post_owf; kwargs...)
end

function post_owf(wm::AbstractWaterModel)
    function_head_loss(wm)

    variable_reservoir(wm)
    variable_tank(wm)
    variable_check_valve(wm)
    variable_head(wm)
    variable_flow(wm)
    variable_volume(wm)
    variable_pump(wm)

    for (a, pipe) in ref(wm, :pipe)
        constraint_link_flow(wm, a)

        # TODO: Call this something other than status.
        if pipe["status"] == "CV"
            constraint_check_valve(wm, a)
            constraint_head_loss_check_valve(wm, a)
        else
            constraint_head_loss_pipe(wm, a)
        end
    end

    for a in ids(wm, :pump)
        constraint_link_flow(wm, a)
        constraint_head_gain_pump(wm, a)
    end

    for (i, node) in ref(wm, :node)
        constraint_flow_conservation(wm, i)

        #if junction["demand"] > 0.0
        #    constraint_sink_flow(wm, i)
        #end
    end

    #for i in collect(ids(wm, :reservoir))
    #    constraint_source_flow(wm, i)
    #end

    for i in ids(wm, :tank)
        constraint_link_volume(wm, i)
        constraint_tank_state(wm, i)
    end

    objective_owf(wm)
end

function run_mn_owf(file, model_constructor, optimizer; kwargs...)
    return run_model(file, model_constructor, optimizer, post_mn_owf; multinetwork=true, kwargs...)
end

function post_mn_owf(wm::AbstractWaterModel)
    function_head_loss(wm)

    for (n, network) in nws(wm)
        variable_reservoir(wm, n)
        variable_tank(wm, n)
        variable_check_valve(wm, n)
        variable_head(wm, n)
        variable_flow(wm, n)
        variable_volume(wm, n)
        variable_pump(wm, n)

        for (a, pipe) in ref(wm, :pipe, nw=n)
            constraint_link_flow(wm, a, nw=n)

            # TODO: Call this something other than status.
            if pipe["status"] == "CV"
                constraint_check_valve(wm, a, nw=n)
                constraint_head_loss_check_valve(wm, a, nw=n)
            else
                constraint_head_loss_pipe(wm, a, nw=n)
            end
        end

        for a in ids(wm, :pump, nw=n)
            constraint_link_flow(wm, a, nw=n)
            constraint_head_gain_pump(wm, a, nw=n)
        end

        for (i, node) in ref(wm, :node, nw=n)
            constraint_flow_conservation(wm, i, nw=n)

            #if junction["demand"] > 0.0
            #    constraint_sink_flow(wm, i, nw=n)
            #end
        end

        #for i in collect(ids(wm, n, :reservoir))
        #    constraint_source_flow(wm, i, nw=n)
        #end
 
        for i in ids(wm, :tank, nw=n)
            constraint_link_volume(wm, i, nw=n)
        end
    end

    network_ids = sort(collect(nw_ids(wm)))

    n_1 = network_ids[1]

    # Initial conditions of tanks.
    for i in ids(wm, :tank, nw=n_1)
        constraint_tank_state(wm, i, nw=n_1)
    end

    # Pump and tank states after the initial time step.
    for n_2 in network_ids[2:end]
        for i in ids(wm, :tank, nw=n_2)
            constraint_tank_state(wm, i, n_1, n_2)
        end

        n_1 = n_2
    end

    objective_owf(wm)
end
