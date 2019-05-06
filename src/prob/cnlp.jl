export run_cnlp_hw

function run_cnlp_hw(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_cnlp_hw; kwargs...)
end

function run_cnlp_hw(file, modifications_path, model_constructor, solver; kwargs...)
    return run_generic_model(file, modifications_path, model_constructor, solver, post_cnlp_hw; kwargs...)
end

function post_cnlp_hw(wm::GenericWaterModel; kwargs...)
    variable_directed_flow(wm, wm.cnw)

    for i in collect(ids(wm, wm.cnw, :junctions))
        constraint_directed_flow_conservation(wm, i, wm.cnw)
    end

    objective_cnlp_hw(wm)
end
