export run_feasibility_relaxed, run_feasibility_exact

function run_feasibility_relaxed(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_feasibility_relaxed; kwargs...)
end

function run_feasibility_exact(file, model_constructor, solver_relaxed, solver_exact; kwargs...)
    # Solve the relaxed problem and save directional data.
    wm_relaxed = build_generic_model(file, model_constructor, post_feasibility_relaxed; kwargs...)
    relaxed_solution = solve_generic_model(wm_relaxed, solver_relaxed; solution_builder = get_solution)
    update_flow_directions(wm_relaxed.data, relaxed_solution["solution"])

    # Solve and return the exact model.
    wm_exact = build_generic_model(wm_relaxed.data, model_constructor, post_feasibility_exact; kwargs...)
    exact_solution = solve_generic_model(wm_exact, solver_exact; solution_builder = get_solution)
    return exact_solution
end

function post_feasibility_relaxed(wm::GenericWaterModel; kwargs...)
    num_separators = wm.setting["num_separators"]

    variable_flow(wm)
    variable_head(wm)
    variable_head_difference(wm)
    variable_flow_direction(wm)

    for i in [collect(ids(wm, :junctions)); collect(ids(wm, :reservoirs))]
        constraint_flow_conservation(wm, i)
    end

    for a in collect(ids(wm, :pipes))
        constraint_potential_flow_coupling(wm, a, true, num_separators)
        constraint_define_gamma(wm, a)
        constraint_flow_direction(wm, a)
    end

    # Uncomment this to get physical feasibility immediately.
    # objective_minimize_gamma(wm)
end

function post_feasibility_exact(wm::GenericWaterModel; kwargs...)
    variable_flow(wm)
    variable_head(wm)
    variable_head_difference(wm)

    for i in [collect(ids(wm, :junctions)); collect(ids(wm, :reservoirs))]
        constraint_flow_conservation(wm, i)
    end

    for a in collect(ids(wm, :pipes))
        constraint_potential_flow_coupling(wm, a, false)
        constraint_define_gamma(wm, a)
    end
end
