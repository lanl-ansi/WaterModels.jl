export run_cvx_hw, run_cvx_dw

function run_cvx_hw(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_cvx_hw; kwargs...)
end

function run_cvx_hw(file, modifications_path, model_constructor, solver; kwargs...)
    return run_generic_model(file, modifications_path, model_constructor, solver, post_cvx_hw; kwargs...)
end

function run_cvx_dw(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_cvx_dw; kwargs...)
end

function run_cvx_dw(file, modifications_path, model_constructor, solver; kwargs...)
    return run_generic_model(file, modifications_path, model_constructor, solver, post_cvx_dw; kwargs...)
end

function post_cvx_hw(wm::GenericWaterModel; kwargs...)
    variable_directed_flow(wm)

    for i in collect(ids(wm, :junctions))
        constraint_directed_flow_conservation(wm, i)
    end

    objective_cvx_hw(wm)
end

function post_cvx_dw(wm::GenericWaterModel; kwargs...)
    #variable_flow_cvxnlp(wm)

    #for i in collect(ids(wm, :junctions))
    #    constraint_flow_conservation_cvx(wm, i)
    #end

    #objective_cvxnlp(wm, 2.0)
end
