export run_wf

function run_wf(network, model_constructor, optimizer; kwargs...)
    post_wf = get_post_wf(model_constructor, kwargs[:alpha])
    return run_generic_model(network, model_constructor, optimizer, post_wf; kwargs...)
end

#function get_post_wf(T::AbstractWaterFormulation, alpha::Float64, n::Int=wm.cnw; kwargs...)
#    if T == StandardMILPRForm
#        return function (wm::GenericWaterModel{T}, n::Int=wm.cnw; kwargs...)
#            variable_head(wm)
#            variable_directed_flow(wm)
#
#            variable_head_difference(wm)
#            variable_flow_direction(wm)
#
#            for a in collect(ids(wm, :links))
#                constraint_flow_direction_selection(wm, a)
#                constraint_head_difference(wm, a)
#                constraint_potential_loss(wm, a)
#                constraint_potential_loss_slope(wm, a)
#            end
#
#            for i in collect(ids(wm, :reservoirs))
#                constraint_source_flow(wm, i)
#            end
#
#            for (i, junction) in wm.ref[:nw][n][:junctions]
#                constraint_directed_flow_conservation(wm, i)
#
#                if junction["demand"] > 0.0
#                    constraint_sink_flow(wm, i)
#                end
#            end
#        end
#    elseif T == StandardCVXNLPForm
#        return function (wm::GenericWaterModel{T}, n::Int=wm.cnw; kwargs...)
#            function_if_alpha(wm, alpha)
#            variable_directed_flow(wm)
#
#            for (i, junction) in wm.ref[:nw][n][:junctions]
#                constraint_directed_flow_conservation(wm, i)
#            end
#        end
#    end
#end

function get_post_wf(alpha::Float64; kwargs...)
    return function (wm::GenericWaterModel{T}, n::Int=wm.cnw; kwargs...) where T <: AbstractWaterFormulation
        if T <: Union{AbstractCVXNLPForm, AbstractMINLPForm, AbstractMICPForm}
            function_if_alpha(wm, alpha)
        end

        variable_flow(wm, n, alpha=alpha)
        variable_head(wm, n)

        for a in collect(ids(wm, :links))
            constraint_potential_loss(wm, a, n)
        end

        for i in collect(ids(wm, :junctions))
            constraint_flow_conservation(wm, i, n)
        end

        objective_wf(wm)

        #variable_head_difference(wm)
        #variable_directed_flow(wm)
        #variable_flow_direction(wm)

        #for a in collect(ids(wm, :links))
        #    constraint_flow_direction_selection(wm, a)
        #    constraint_head_difference(wm, a)
        #    constraint_potential_loss(wm, a)
        #    constraint_potential_loss_slope(wm, a)
        #end

        #for i in collect(ids(wm, :reservoirs))
        #    constraint_source_flow(wm, i)
        #end

        #for (i, junction) in wm.ref[:nw][n][:junctions]
        #    constraint_directed_flow_conservation(wm, i)

        #    if junction["demand"] > 0.0
        #        constraint_sink_flow(wm, i)
        #    end
        #end
    end
end
