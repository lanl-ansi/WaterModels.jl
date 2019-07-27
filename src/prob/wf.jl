function run_wf(network, model_constructor, optimizer; relaxed::Bool=false, kwargs...)
    return run_generic_model(network, model_constructor, optimizer, post_wf, relaxed=relaxed; kwargs...)
end

function post_wf(wm::GenericWaterModel{T}) where T
    # TODO: Do not use conditionals here.
    if T <: Union{AbstractMICPForm, AbstractNCNLPForm}
        function_f_alpha(wm, convex=false)
    elseif T <: AbstractCNLPForm
        function_if_alpha(wm, convex=true)
    end

    variable_reservoir(wm)
    variable_tank(wm)
    variable_check_valve(wm)
    variable_head(wm)
    variable_flow(wm)
    variable_volume(wm)
    variable_pump(wm)

    for (a, pipe) in ref(wm, :pipes)
        constraint_link_flow(wm, a)

        # TODO: Call this something other than status.
        if pipe["status"] == "CV"
            constraint_check_valve(wm, a)
            constraint_potential_loss_check_valve(wm, a)
        else
            constraint_potential_loss_pipe(wm, a)
        end
    end

    for a in ids(wm, :pumps)
        constraint_link_flow(wm, a)
        constraint_potential_loss_pump(wm, a)
        constraint_pump_control(wm, a, initial=true)
    end

    for (i, node) in ref(wm, :nodes)
        constraint_flow_conservation(wm, i)

        #if junction["demand"] > 0.0
        #    constraint_sink_flow(wm, i)
        #end
    end

    ##for i in collect(ids(wm, :reservoirs))
    ##    constraint_source_flow(wm, i)
    ##end

    for i in ids(wm, :tanks)
        constraint_link_volume(wm, i)
        constraint_tank_state(wm, i)
    end

    objective_wf(wm)
end

function run_mn_wf(file, model_constructor, optimizer; kwargs...)
    return run_generic_model(file, model_constructor, optimizer, post_mn_wf; multinetwork=true, kwargs...)
end

function post_mn_wf(wm::GenericWaterModel{T}) where T
    if T <: Union{AbstractMICPForm, AbstractNCNLPForm}
        function_f_alpha(wm, convex=false)
    elseif T <: AbstractCNLPForm
        function_if_alpha(wm, convex=true)
    end

    for (n, network) in nws(wm)
        variable_reservoir(wm, n)
        variable_tank(wm, n)
        variable_check_valve(wm, n)
        variable_head(wm, n)
        variable_flow(wm, n)
        variable_volume(wm, n)
        variable_pump(wm, n)

        for (a, pipe) in ref(wm, :pipes)
            constraint_link_flow(wm, a, nw=n)

            # TODO: Call this something other than status.
            if pipe["status"] == "CV"
                constraint_check_valve(wm, a, nw=n)
                constraint_potential_loss_check_valve(wm, a, nw=n)
            else
                constraint_potential_loss_pipe(wm, a, nw=n)
            end
        end

        for a in ids(wm, :pumps, nw=n)
            constraint_link_flow(wm, a, nw=n)
            constraint_potential_loss_pump(wm, a, nw=n)
        end

        for (i, node) in ref(wm, :nodes, nw=n)
            constraint_flow_conservation(wm, i, nw=n)

            #if junction["demand"] > 0.0
            #    constraint_sink_flow(wm, i, nw=n)
            #end
        end

        #for i in collect(ids(wm, n, :reservoirs))
        #    constraint_source_flow(wm, i, nw=n)
        #end
 
        for i in ids(wm, n, :tanks)
            constraint_link_volume(wm, i, nw=n)
        end
    end

    network_ids = sort(collect(nw_ids(wm)))

    n_1 = network_ids[1]

    # Initial conditions of pumps.
    for a in ids(wm, :pumps, nw=n_1)
        constraint_pump_control(wm, a, nw=n_1, initial=true)
    end

    # Initial conditions of tanks.
    for i in ids(wm, :tanks, nw=n_1)
        constraint_tank_state(wm, i, nw=n_1)
    end

    # Pump and tank states after the initial time step.
    for n_2 in network_ids[2:end]
       for i in ids(wm, :tanks, nw=n_2)
           constraint_tank_state(wm, i, n_1, n_2)
       end

        for a in ids(wm, :pumps, nw=n_2)
            constraint_pump_control(wm, a, n_1, n_2)
        end

       n_1 = n_2
    end

    objective_wf(wm)
end
