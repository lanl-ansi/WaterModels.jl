export run_ne

function run_ne(network, model_constructor, optimizer; relaxed::Bool=false, kwargs...)
    return run_generic_model(network, model_constructor, optimizer, post_ne, relaxed=relaxed; kwargs...)
end

function post_ne(wm::GenericWaterModel{T}, n::Int=wm.cnw; kwargs...) where T <: AbstractWaterFormulation
    if T <: AbstractNCNLPForm
        function_f_alpha(wm, n, convex=false)
    elseif T <: AbstractMICPForm
        function_f_alpha(wm, n, convex=true)
    elseif T <: AbstractCNLPForm
        Memento.error(_LOGGER, "CNLP formulation does not support network expansion.")
    end

    variable_head(wm, n)
    variable_flow(wm, n)
    variable_flow_ne(wm, n)
    variable_resistance_ne(wm, n)

    for a in ids(wm, n, :links)
        constraint_link_flow(wm, a, n)
    end

    for a in setdiff(ids(wm, n, :links), ids(wm, n, :links_ne))
        constraint_potential_loss(wm, a, n)
    end

    for a in collect(ids(wm, n, :links_ne))
        constraint_potential_loss_ne(wm, a, n)
        constraint_resistance_selection_ne(wm, a, n)
        constraint_link_flow_ne(wm, a, n)
    end

    for (i, junction) in ref(wm, n, :junctions)
        constraint_flow_conservation(wm, i, n)

        if junction["demand"] > 0.0
            constraint_sink_flow(wm, i, n)
        end
    end

    for i in collect(ids(wm, n, :reservoirs))
        constraint_source_flow(wm, i, n)
    end

    objective_ne(wm, n)
end
