export run_ne_hw, run_ne_dw, post_ne_hw, post_ne_hw_segmented, post_ne_dw

function run_ne_hw(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_ne_hw; kwargs...)
end

function run_ne_hw(file, modifications_path, model_constructor, solver; kwargs...)
    return run_generic_model(file, modifications_path, model_constructor, solver, post_ne_hw; kwargs...)
end

function run_ne_dw(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_ne_dw; kwargs...)
end

function run_ne_dw(file, modifications_path, model_constructor, solver; kwargs...)
    return run_generic_model(file, modifications_path, model_constructor, solver, post_ne_dw; kwargs...)
end

function post_ne_hw_segmented(wm::GenericWaterModel, n_n::Int = wm.cnw, n_s::Int = 2)
    variable_head(wm, n_n)
    variable_segment(wm, n_n, n_s)
    variable_segmented_directed_flow(wm, n_n, n_s)

    variable_head_difference(wm, n_n)
    variable_flow_direction(wm, n_n)
    variable_resistance(wm, n_n)

    for a in collect(ids(wm, :connection))
        constraint_select_segment(wm, a, n_n)
        constraint_select_resistance(wm, a, n_n)
        constraint_select_segmented_flow_term(wm, a, n_n, n_s)
        constraint_head_difference(wm, a, n_n)
        #constraint_potential_loss_segmented(wm, a, n_n, n_s)
        #constraint_potential_loss_slope_segmented(wm, a, n_n, n_s)
    end

    for (i, junction) in wm.ref[:nw][n_n][:junctions]
        constraint_directed_segmented_flow_conservation(wm, i, n_n)
    end

    objective_minimize_resistance_cost(wm)
end

function post_ne_hw(wm::GenericWaterModel, n_n::Int = wm.cnw; kwargs...)
    variable_head(wm, n_n)
    variable_directed_flow(wm, n_n)

    variable_head_difference(wm, n_n)
    variable_flow_direction(wm, n_n)
    variable_resistance(wm, n_n)

    for a in collect(ids(wm, :connection))
        constraint_select_resistance(wm, a, n_n)
        constraint_select_flow_term(wm, a, n_n)
        constraint_head_difference(wm, a, n_n)
        constraint_potential_loss(wm, a, n_n)
        constraint_potential_loss_slope(wm, a, n_n)
    end

    for (i, junction) in wm.ref[:nw][wm.cnw][:junctions]
        constraint_directed_flow_conservation(wm, i, n_n)

        #if junction["demand"] > 0.0
        #    constraint_sink_flow(wm, i)
        #end
    end

    #for i in collect(ids(wm, :reservoirs))
    #    constraint_source_flow(wm, i)
    #end

    objective_minimize_resistance_cost(wm)
end

function post_ne_bt_hw(wm::GenericWaterModel; kwargs...)
    variable_flow(wm)
    variable_head_ne(wm)
    variable_pipe_ne(wm)

    for i in [collect(ids(wm, :junctions)); collect(ids(wm, :reservoirs))]
        constraint_junction_mass_flow(wm, i)
    end

    for a in collect(ids(wm, :connection_unknown_direction))
        constraint_hw_unknown_direction_ne(wm, a)
    end

    for a in collect(ids(wm, :connection_known_direction))
        constraint_hw_known_direction(wm, a)
    end

    #variable_objective_ne(wm)
    #h_i = wm.var[:nw][wm.cnw][:h]["2"]
    #objective_minimize_variable(wm, h_i)
end

function post_ne_dw(wm::GenericWaterModel; kwargs...)
    variable_flow(wm)
    variable_head(wm)

    for i in [collect(ids(wm, :junctions)); collect(ids(wm, :reservoirs))]
        constraint_flow_conservation(wm, i)
    end

    for a in collect(ids(wm, :connection_unknown_direction))
        constraint_dw_unknown_direction(wm, a)
    end

    for a in collect(ids(wm, :connection_known_direction))
        constraint_dw_known_direction(wm, a)
    end
end
