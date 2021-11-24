function solve_des(network, model_constructor, optimizer; kwargs...)
    return solve_model(network, model_constructor, optimizer, build_des; kwargs...)
end


function build_des(wm::AbstractWaterModel)
    # Build the water flow problem.
    build_wf(wm)

    # Add the network design objective.
    objective_des(wm)
end
