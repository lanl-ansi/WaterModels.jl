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

    for i in [collect(ids(wm, :junctions)); collect(ids(wm, :reservoirs))]
        constraint_flow_conservation(wm, i)
        #if length(wm.ref[:nw][wm.cnw][:junction_connections][i]) == 2 &&
        #    haskey(wm.ref[:nw][wm.cnw][:junctions], i)
        #    if wm.ref[:nw][wm.cnw][:junctions][i]["demand"] == 0.0
        #        constraint_degree_two(wm, i)
        #    end
        #end
    end

    for a in collect(ids(wm, :pipes))
        constraint_hw_unknown_direction(wm, a)
    end
end

function post_wf_dw(wm::GenericWaterModel; kwargs...)
    variable_flow(wm)
    variable_head(wm)

    for i in [collect(ids(wm, :junctions)); collect(ids(wm, :reservoirs))]
        constraint_flow_conservation(wm, i)
        #if length(wm.ref[:nw][wm.cnw][:junction_connections][i]) == 2 &&
        #    haskey(wm.ref[:nw][wm.cnw][:junctions], i)
        #    if wm.ref[:nw][wm.cnw][:junctions][i]["demand"] == 0.0
        #        constraint_degree_two(wm, i)
        #    end
        #end
    end

    for a in collect(ids(wm, :pipes))
        constraint_dw_unknown_direction(wm, a)
    end
end

#export run_wf_relaxed, run_wf_exact
#
#" entry point into running the gas flow feasability problem "
#function run_wf(file, model_constructor, solver; kwargs...)
#    return run_generic_model(file, model_constructor, solver, post_wf; kwargs...) 
#end
#
#function run_wf_relaxed(file, model_constructor, solver; kwargs...)
#    return run_generic_model(file, model_constructor, solver, post_wf_relaxed; kwargs...)
#end

#function run_wf_exact(file, model_constructor, solver_relaxed, solver_exact; kwargs...)
#    exact_solution = Dict{String,Any}()
#    status = :LocalInfeasibile
#
#    # Build the relaxed model.
#    wm_relaxed = build_generic_model(file, model_constructor, post_wf_relaxed; kwargs...)
#
#    while status == :LocalInfeasibile || status == :Infeasible
#        # Solve the relaxed model.
#        setsolver(wm_relaxed.model, solver_relaxed)
#        status, solve_time = solve(wm_relaxed)
#        update_flow_directions(wm_relaxed.data, wm_relaxed) #relaxed_solution["solution"])
#
#        # Solve and return the exact model.
#        wm_exact = build_generic_model(wm_relaxed.data, model_constructor, post_wf_exact; kwargs...)
#        setsolver(wm_exact.model, solver_exact)
#        status, solve_time = solve(wm_exact)
#
#        # If the solution is infeasible, the directions are infeasible.
#        if status == :LocalInfeasible || status == :Infeasible
#            constraint_no_good(wm_relaxed)
#            reset_flow_directions(wm_relaxed.data)
#        else
#            exact_solution = build_solution(wm_exact, status, solve_time; solution_builder = get_solution)
#        end
#    end
#
#    return exact_solution
#end
#
#function post_wf_relaxed(wm::GenericWaterModel; kwargs...)
#    num_separators = wm.setting["num_separators"]
#
#    variable_flow(wm)
#    variable_head(wm)
#    variable_head_difference(wm)
#    variable_flow_direction(wm)
#
#    for i in [collect(ids(wm, :junctions)); collect(ids(wm, :reservoirs))]
#        constraint_flow_conservation(wm, i)
#
#        if length(wm.ref[:nw][wm.cnw][:junction_connections][i]) == 2 &&
#            haskey(wm.ref[:nw][wm.cnw][:junctions], i)
#            if wm.ref[:nw][wm.cnw][:junctions][i]["demand"] == 0.0
#                constraint_degree_two(wm, i)
#            end
#        end
#    end
#
#    for a in collect(ids(wm, :pipes))
#        constraint_potential_flow_coupling(wm, a, true, num_separators)
#        constraint_define_gamma(wm, a)
#        constraint_flow_direction(wm, a)
#    end
#
#    # Uncomment this to get physical feasibility immediately.
#    # objective_minimize_gamma(wm)
#end
#
#function post_wf_exact(wm::GenericWaterModel; kwargs...)
#    variable_flow(wm)
#    variable_head(wm)
#    variable_head_difference(wm)
#
#    for i in [collect(ids(wm, :junctions)); collect(ids(wm, :reservoirs))]
#        constraint_flow_conservation(wm, i)
#    end
#
#    for a in collect(ids(wm, :pipes))
#        constraint_potential_flow_coupling(wm, a, false)
#        constraint_define_gamma(wm, a)
#    end
#end

