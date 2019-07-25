function run_owf(network, model_constructor, optimizer; relaxed::Bool=false, kwargs...)
    return run_generic_model(network, model_constructor, optimizer, post_owf, relaxed=relaxed; kwargs...)
end

function post_owf(wm::GenericWaterModel{T}) where T
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

    for a in setdiff(ids(wm, :pipes), ids(wm, :check_valves))
        constraint_potential_loss_pipe(wm, a)
        constraint_link_flow(wm, a)
    end

    for a in ids(wm, :check_valves)
        constraint_check_valve(wm, a)
        constraint_potential_loss_check_valve(wm, a)
        constraint_link_flow(wm, a)
    end

    for a in ids(wm, :pumps)
        constraint_potential_loss_pump(wm, a)
    end

    for a in ids(wm, :check_valves)
        constraint_check_valve(wm, a)
    end

    for (i, node) in ref(wm, :nodes)
        constraint_flow_conservation(wm, i)

        #if junction["demand"] > 0.0
        #    constraint_sink_flow(wm, i)
        #end
    end

    #for i in collect(ids(wm, :reservoirs))
    #    constraint_source_flow(wm, i)
    #end
    
    for i in ids(wm, :tanks)
        constraint_link_volume(wm, i)
    end

    objective_owf(wm)
end

function run_mn_owf(file, model_constructor, optimizer; relaxed::Bool=false, kwargs...)
    return run_generic_model(file, model_constructor, optimizer, post_mn_owf, relaxed=relaxed; multinetwork=true, kwargs...)
end

function post_mn_owf(wm::GenericWaterModel{T}) where T
    if T <: Union{AbstractMICPForm, AbstractNCNLPForm}
        function_f_alpha(wm, wm.cnw, convex=false)
    elseif T <: AbstractCNLPForm
        function_if_alpha(wm, wm.cnw, convex=true)
    end

    for (n, network) in nws(wm)
        variable_reservoir(wm, n)
        variable_tank(wm, n)
        variable_check_valve(wm, n)
        variable_head(wm, n)
        variable_flow(wm, n)
        variable_volume(wm, n)
        variable_pump(wm, n)

        for a in setdiff(ids(wm, n, :pipes), ids(wm, n, :check_valves))
            constraint_potential_loss_pipe(wm, a, n)
            constraint_link_flow(wm, a, nw=n)
        end

        for a in ids(wm, n, :check_valves)
            constraint_check_valve(wm, a, n)
            constraint_potential_loss_check_valve(wm, a, n)
            constraint_link_flow(wm, a, n)
        end

        for a in ids(wm, n, :pumps)
            constraint_potential_loss_pump(wm, a, n)
        end

        for (i, node) in  ref(wm, n, :nodes)
            constraint_flow_conservation(wm, i, nw=n)

            #if junction["demand"] > 0.0
            #    constraint_sink_flow(wm, i, nw=n)
            #end
        end

        #for i in collect(ids(wm, n, :reservoirs))
        #    constraint_source_flow(wm, i, nw=n)
        #end
        
        for i in ids(wm, n, :tanks)
            constraint_link_volume(wm, i, n)
        end
    end

    network_ids = sort(collect(nw_ids(wm)))

    n_1 = network_ids[1]

    for i in ids(wm, :tanks, nw=n_1)
        constraint_tank_state(wm, i, nw=n_1)
    end

    for n_2 in network_ids[2:end]
       for i in ids(wm, :tanks, nw=n_2)
           constraint_tank_state(wm, i, n_1, n_2)
       end

       n_1 = n_2
    end

    objective_owf(wm)
end
