export run_wf_hw

function run_wf_hw(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_wf_hw; kwargs...)
end

function run_wf_hw(file, modifications_path, model_constructor, solver; kwargs...)
    return run_generic_model(file, modifications_path, model_constructor, solver, post_wf_hw; kwargs...)
end

function post_wf_hw(wm::GenericWaterModel; kwargs...)
    variable_head(wm)
    variable_directed_flow(wm)

    variable_head_difference(wm)
    variable_flow_direction(wm)
    variable_resistance(wm)

    for a in collect(ids(wm, :connection))
        constraint_select_resistance(wm, a)
        constraint_select_flow_term(wm, a)
        constraint_head_difference(wm, a)
        constraint_potential_loss(wm, a)
        constraint_potential_loss_slope(wm, a)
    end

    for i in collect(ids(wm, :junctions))
        constraint_directed_flow_conservation(wm, i)
    end
end
