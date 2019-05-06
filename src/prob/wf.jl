export run_wf

function run_wf(network, model_constructor, optimizer; relaxed::Bool=false, kwargs...)
    post_wf = get_post_wf(kwargs[:alpha])
    return run_generic_model(network, model_constructor, optimizer, post_wf, relaxed=relaxed; kwargs...)
end

function get_post_wf(alpha::Float64; kwargs...)
    function (wm::GenericWaterModel{T}, n::Int=wm.cnw; kwargs...) where T <: AbstractWaterFormulation
        if T <: Union{AbstractMICPForm, AbstractNCNLPForm}
            function_f_alpha(wm, alpha - 1.0)
        elseif T <: AbstractCNLPForm
            function_if_alpha(wm, alpha - 1.0)
        end

        variable_head(wm, n)
        variable_flow(wm, n, alpha=alpha)

        for a in collect(ids(wm, n, :links))
            constraint_potential_loss(wm, a, n, alpha=alpha)
            constraint_link_flow(wm, a, n)
        end

        for (i, junction) in wm.ref[:nw][n][:junctions]
            constraint_flow_conservation(wm, i, n)

            if junction["demand"] > 0.0
                constraint_sink_flow(wm, i, n)
            end
        end

        for i in collect(ids(wm, n, :reservoirs))
            constraint_source_flow(wm, i, n)
        end

        objective_wf(wm, n)
    end
end
