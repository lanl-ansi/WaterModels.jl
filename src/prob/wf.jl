export run_wf

function run_wf(file, model_constructor, optimizer; kwargs...)
    post_wf = get_post_wf(kwargs[:alpha]; kwargs...)
    return run_generic_model(file, model_constructor, optimizer, post_wf; kwargs...)
end

function get_post_wf(alpha::Float64; kwargs...)
    return function (wm::GenericWaterModel{T}; kwargs...) where T <: StandardCVXNLPForm
        function_f_alpha(wm, alpha)
        function_if_alpha(wm, alpha)

        variable_directed_flow(wm)
        variable_wf_objective_term(wm)

        for a in collect(ids(wm, :connection))
            constraint_wf_objective_term(wm, a)
        end

        for i in collect(ids(wm, :junctions))
            constraint_directed_flow_conservation(wm, i)
        end

        objective_wf(wm)
    end
end
#
#function post_wf(wm::GenericWaterModel; kwargs...)
#    function_f_alpha(wm, alpha = alpha)
#    function_if_alpha(wm, alpha = alpha)
#
#
#    #for i in collect(ids(wm, :junctions))
#    #    constraint_directed_flow_conservation(wm, i)
#    #end
#
#    #objective_cvx(wm)
#end

#function post_wf_hw(wm::GenericWaterModel; kwargs...)
#    variable_head(wm)
#    variable_directed_flow(wm)
#
#    variable_head_difference(wm)
#    variable_flow_direction(wm)
#    variable_resistance(wm)
#
#    for a in collect(ids(wm, :connection))
#        constraint_select_resistance(wm, a)
#        constraint_select_flow_term(wm, a)
#        constraint_head_difference(wm, a)
#        constraint_potential_loss(wm, a)
#        constraint_potential_loss_slope(wm, a)
#    end
#
#    for i in collect(ids(wm, :junctions))
#        constraint_directed_flow_conservation(wm, i)
#    end
#end
#
#function post_wf_dw(wm::GenericWaterModel; kwargs...)
#    #variable_flow(wm)
#    #variable_head(wm)
#
#    #for i in [collect(ids(wm, :junctions)); collect(ids(wm, :reservoirs))]
#    #    constraint_junction_mass_flow(wm, i)
#    #end
#
#    #for a in collect(ids(wm, :connection_unknown_direction))
#    #    constraint_dw_unknown_direction(wm, a)
#    #end
#
#    #for a in collect(ids(wm, :connection_known_direction))
#    #    constraint_dw_known_direction(wm, a)
#    #end
#
#    #objective_dummy(wm)
#end
