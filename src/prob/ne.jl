export run_ne_hw, run_ne_dw

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

function post_ne_hw(wm::GenericWaterModel; kwargs...)
    variable_flow(wm)
    variable_head(wm)
    variable_pipe_ne(wm)

    for i in [collect(ids(wm, :junctions)); collect(ids(wm, :reservoirs))]
        constraint_flow_conservation(wm, i)
    end

    for a in collect(ids(wm, :connection_unknown_direction))
        constraint_hw_unknown_direction_ne(wm, a)
    end

    for a in collect(ids(wm, :connection_known_direction))
        constraint_hw_known_direction(wm, a)
    end
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
