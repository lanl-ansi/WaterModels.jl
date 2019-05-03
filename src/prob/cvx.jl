export run_cvx_hw

function run_cvx_hw(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_cvx_hw; kwargs...)
end

function run_cvx_hw(file, modifications_path, model_constructor, solver; kwargs...)
    return run_generic_model(file, modifications_path, model_constructor, solver, post_cvx_hw; kwargs...)
end

function post_cvx_hw(wm::GenericWaterModel; kwargs...)
    variable_directed_flow(wm, wm.cnw)

    for i in collect(ids(wm, wm.cnw, :junctions))
        constraint_directed_flow_conservation(wm, i, wm.cnw)
    end

    objective_cvx_hw(wm)
end
