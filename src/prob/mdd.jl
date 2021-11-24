function solve_mdd(network, model_constructor, optimizer; kwargs...)
    return solve_model(network, model_constructor, optimizer, build_mdd; kwargs...)
end


function solve_mn_mdd(file, model_constructor, optimizer; kwargs...)
    return solve_model(file, model_constructor, optimizer, build_mn_mdd; multinetwork=true, kwargs...)
end


function build_mdd(wm::AbstractWaterModel)
    # Build the water flow problem.
    build_wf(wm)

    # Add the demand maximization objective.
    objective_max_demand(wm)
end


function build_mn_mdd(wm::AbstractWaterModel)
    # Build the water flow problem.
    build_mn_wf(wm)

    # Add the demand maximization objective.
    objective_max_demand(wm)
end