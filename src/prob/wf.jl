export run_wf_hw, run_wf_dw

function run_wf_hw(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_wf_hw; kwargs...)
end

function run_wf_hw(file, modifications_path, model_constructor, solver; kwargs...)
    return run_generic_model(file, modifications_path, model_constructor, solver, post_wf_hw; kwargs...)
end

function run_wf_dw(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_wf_dw; kwargs...)
end

function run_wf_dw(file, modifications_path, model_constructor, solver; kwargs...)
    return run_generic_model(file, modifications_path, model_constructor, solver, post_wf_dw; kwargs...)
end

function post_wf_hw(wm::GenericWaterModel; kwargs...)
    variable_head(wm)
    variable_directed_flow(wm)

    variable_head_difference(wm)
    variable_flow_direction(wm)
    variable_resistance(wm)

    function_head_loss_hw(wm)

    for a in collect(ids(wm, :connection))
        constraint_select_resistance(wm, a)
        constraint_select_flow_term(wm, a)
        constraint_head_difference(wm, a)
        constraint_potential_loss(wm, a)
        #constraint_potential_loss_slope(wm, a)
    end

    for i in collect(ids(wm, :junctions))
        constraint_directed_flow_conservation(wm, i)
    end

    objective_minimize_resistance_cost(wm)
end

function post_wf_dw(wm::GenericWaterModel; kwargs...)
    #variable_flow(wm)
    #variable_head(wm)

    #for i in [collect(ids(wm, :junctions)); collect(ids(wm, :reservoirs))]
    #    constraint_junction_mass_flow(wm, i)
    #end

    #for a in collect(ids(wm, :connection_unknown_direction))
    #    constraint_dw_unknown_direction(wm, a)
    #end

    #for a in collect(ids(wm, :connection_known_direction))
    #    constraint_dw_known_direction(wm, a)
    #end

    #objective_dummy(wm)
end
