function solve_mdd(network, model_constructor, optimizer; kwargs...)
    return solve_model(network, model_constructor, optimizer, build_mdd; kwargs...)
end


function run_mdd(network, model_constructor, optimizer; kwargs...)
    Memento.warn(_LOGGER, "\"run_\" methods should be renamed \"solve_\" and will be deprecated in future versions.")
    return solve_owf(network, model_constructor, optimizer; kwargs...)
end


function solve_mn_mdd(file, model_constructor, optimizer; kwargs...)
    return solve_model(file, model_constructor, optimizer, build_mn_mdd; multinetwork=true, kwargs...)
end


function run_mn_mdd(network, model_constructor, optimizer; kwargs...)
    Memento.warn(_LOGGER, "\"run_\" methods should be renamed \"solve_\" and will be deprecated in future versions.")
    return solve_mn_mdd(network, model_constructor, optimizer; kwargs...)
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