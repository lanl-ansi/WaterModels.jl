export run_feasibility

""
function run_feasibility(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_feasibility; kwargs...)
end

""
function post_feasibility(wm::GenericWaterModel)
    variable_flow(wm)
end
