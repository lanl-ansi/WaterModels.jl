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
    variable_head(wm)
    variable_flow(wm)
    variable_pump(wm)

    for a in collect(ids(wm, :pipes))
        constraint_potential_loss_pipe(wm, a)
        constraint_link_flow(wm, a)
    end

    for a in collect(ids(wm, :pumps))
        constraint_potential_loss_pump(wm, a)
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

    objective_wf(wm)
end

function run_mn_wf(file, model_constructor, optimizer; kwargs...)
    return run_generic_model(file, model_constructor, optimizer, post_mn_wf; multinetwork=true, kwargs...)
end

function post_mn_wf(wm::GenericWaterModel{T}) where T
    for (n, network) in nws(wm)
        if T <: Union{AbstractMICPForm, AbstractNCNLPForm}
            function_f_alpha(wm, n, convex=false)
        elseif T <: AbstractCNLPForm
            function_if_alpha(wm, n, convex=true)
        end

        variable_reservoir(wm, n)
        variable_head(wm, n)
        variable_flow(wm, n)
        variable_pump(wm, n)

        for a in collect(ids(wm, n, :pipes))
            constraint_potential_loss_pipe(wm, a, n)
            constraint_link_flow(wm, a, n)
        end

        for a in collect(ids(wm, n, :pumps))
            constraint_potential_loss_pump(wm, a, n)
        end

        for (i, node) in ref(wm, n, :nodes)
            constraint_flow_conservation(wm, i, n)

            #if junction["demand"] > 0.0
            #    constraint_sink_flow(wm, i, n)
            #end
        end

        #for i in collect(ids(wm, n, :reservoirs))
        #    constraint_source_flow(wm, i, n)
        #end
    end

    objective_wf(wm)
end
