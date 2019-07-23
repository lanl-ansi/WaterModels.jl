export run_owf, run_mn_owf

function run_owf(network, model_constructor, optimizer; relaxed::Bool=false, kwargs...)
    return run_generic_model(network, model_constructor, optimizer, post_owf, relaxed=relaxed; kwargs...)
end

function post_owf(wm::GenericWaterModel{T}, n::Int=wm.cnw; kwargs...) where T
    if T <: Union{AbstractMICPForm, AbstractNCNLPForm}
        function_f_alpha(wm, n, convex=false)
    elseif T <: AbstractCNLPForm
        function_if_alpha(wm, n, convex=true)
    end

    variable_head(wm, n)
    variable_flow(wm, n)
    variable_pump(wm, n)

    # TODO: Need something separate for pumps or handle in constraint templates.
    for a in collect(ids(wm, n, :pipes))
        constraint_potential_loss_pipe(wm, a, n)
        constraint_link_flow(wm, a, n)
    end

    for a in collect(ids(wm, n, :pumps))
        constraint_potential_loss_pump(wm, a, n)
    end

    for (i, junction) in ref(wm, n, :junctions)
        constraint_flow_conservation(wm, i, n)

        #if junction["demand"] > 0.0
        #    constraint_sink_flow(wm, i, n)
        #end
    end

    #for i in collect(ids(wm, n, :reservoirs))
    #    constraint_source_flow(wm, i, n)
    #end

    objective_owf(wm)
end

function run_mn_owf(file, model_constructor, optimizer; kwargs...)
    return run_generic_model(file, model_constructor, optimizer, post_mn_owf; multinetwork=true, kwargs...)
end

function post_mn_owf(wm::GenericWaterModel{T}; kwargs...) where T
    if T <: Union{AbstractMICPForm, AbstractNCNLPForm}
        function_f_alpha(wm, wm.cnw, convex=false)
    elseif T <: AbstractCNLPForm
        function_if_alpha(wm, wm.cnw, convex=true)
    end

    for (n, network) in nws(wm)
        variable_head(wm, n)
        variable_flow(wm, n)
        variable_pump(wm, n)

        for a in ids(wm, n, :pipes)
            constraint_potential_loss_pipe(wm, a, n)
            constraint_link_flow(wm, a, n)
        end

        for a in collect(ids(wm, n, :pumps))
            constraint_potential_loss_pump(wm, a, n)
        end

        for (i, junction) in wm.ref[:nw][n][:junctions]
            constraint_flow_conservation(wm, i, n)

            #if junction["demand"] > 0.0
            #    constraint_sink_flow(wm, i, n)
            #end
        end

        #for i in collect(ids(wm, n, :reservoirs))
        #    constraint_source_flow(wm, i, n)
        #end
    end

    objective_owf(wm)
end
