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
    variable_flow(wm)
    variable_head(wm)

    #for i in [collect(ids(wm, :junctions)); collect(ids(wm, :reservoirs))]
    #    constraint_junction_mass_flow(wm, i)
    #end

    #for a in collect(ids(wm, :connection_unknown_direction))
    #    constraint_hw_unknown_direction(wm, a)
    #end

    #for a in collect(ids(wm, :connection_known_direction))
    #    constraint_hw_known_direction(wm, a)
    #end

    objective_dummy(wm)
end

function post_wf_dw(wm::GenericWaterModel; kwargs...)
    variable_flow(wm)
    variable_head(wm)

    for i in [collect(ids(wm, :junctions)); collect(ids(wm, :reservoirs))]
        constraint_junction_mass_flow(wm, i)
    end

    for a in collect(ids(wm, :connection_unknown_direction))
        constraint_dw_unknown_direction(wm, a)
    end

    for a in collect(ids(wm, :connection_known_direction))
        constraint_dw_known_direction(wm, a)
    end

    objective_dummy(wm)
end
