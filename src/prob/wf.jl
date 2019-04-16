export run_wf

function run_wf(network, model_constructor, optimizer; kwargs...)
    post_wf = get_post_wf(kwargs[:alpha]; kwargs...)
    return run_generic_model(network, model_constructor, optimizer, post_wf; kwargs...)
end

function get_post_wf(alpha::Float64; kwargs...)
    return function (wm::GenericWaterModel{T}, n::Int=wm.cnw; kwargs...) where T <: AbstractWaterFormulation
        variable_head(wm)
        variable_directed_flow(wm)

        variable_head_difference(wm)
        variable_flow_direction(wm)

        for a in collect(ids(wm, :links))
            constraint_flow_direction_selection(wm, a)
            constraint_head_difference(wm, a)
            constraint_potential_loss(wm, a)
            constraint_potential_loss_slope(wm, a)
        end

        for i in collect(ids(wm, :reservoirs))
            constraint_source_flow(wm, i)
        end

        for (i, junction) in wm.ref[:nw][n][:junctions]
            constraint_directed_flow_conservation(wm, i)

            if junction["demand"] > 0.0
                constraint_sink_flow(wm, i)
            end
        end
    end
end
