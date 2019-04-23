export run_ne

function run_ne(network, model_constructor, optimizer; kwargs...)
    post_ne = get_post_ne(kwargs[:alpha]; kwargs...)
    return run_generic_model(network, model_constructor, optimizer, post_ne; kwargs...)
end

function get_post_ne(alpha::Float64; kwargs...)
    function (wm::GenericWaterModel{T}, n::Int=wm.cnw; kwargs...) where T <: AbstractWaterFormulation
        if T <: AbstractMINLPForm
            function_f_alpha(wm, alpha)
        end

        variable_head(wm, n)
        variable_flow(wm, n, alpha=alpha)
        variable_resistance_ne(wm, n)

        for (a, arc) in wm.ref[:nw][n][:links]
            constraint_resistance_selection(wm, a, n)
        #    constraint_flow_direction_selection(wm, a, n)
        #    constraint_head_difference(wm, a, n)
        #    constraint_potential_loss(wm, a, n)
        #    constraint_potential_loss_slope(wm, a, n)
        end

        for i in collect(ids(wm, :junctions))
            constraint_flow_conservation(wm, i, n)
        end

        #variable_undirected_flow(wm, n, alpha=alpha)
        #variable_directed_flow(wm, n, alpha=alpha)
        #variable_directed_flow_ne(wm, n, alpha=alpha)
        #variable_directed_flow(wm, n)

        #variable_head_difference(wm, n)
        #variable_flow_direction(wm, n)
        #variable_resistance(wm, n)

        #for (a, arc) in wm.ref[:nw][n][:links]
        #    constraint_resistance_selection(wm, a, n)
        #    constraint_flow_direction_selection(wm, a, n)
        #    constraint_head_difference(wm, a, n)
        #    constraint_potential_loss(wm, a, n)
        #    constraint_potential_loss_slope(wm, a, n)
        #end

        #for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
        #    constraint_source_flow(wm, i)
        #end

        #for (i, junction) in wm.ref[:nw][n][:junctions]
        #    constraint_directed_flow_conservation(wm, i, n)

        #    if junction["demand"] > 0.0
        #        constraint_sink_flow(wm, i)
        #    end
        #end

        #objective_minimize_resistance_cost(wm)
    end
end
