### Optimization-based Bound Tightening for Optimal Water Flow Problems ###
### This is built upon https://github.com/lanl-ansi/PowerModels.jl/blob/master/src/util/obbt.jl


function _fix_indicator(v::JuMP.VariableRef, value::Float64)
    if JuMP.is_binary(v)
        JuMP.fix(v, value)
    else
        JuMP.fix(v, value; force = true)
    end
end


function _unfix_indicator(v::JuMP.VariableRef)
    if JuMP.is_binary(v)
        JuMP.unfix(v)
    else
        JuMP.unfix(v)
        JuMP.set_lower_bound(v, 0.0)
        JuMP.set_upper_bound(v, 1.0)
    end
end


function _solve_bound_problem!(wm::AbstractWaterModel, bound_problem::BoundProblem)
    v = _get_variable_from_index(wm, bound_problem.variable_to_tighten)
    v_one = _get_variable_from_index.(Ref(wm), bound_problem.variables_fix_one)
    v_zero = _get_variable_from_index.(Ref(wm), bound_problem.variables_fix_zero)

    # Determine if the variable being tightened is discrete.
    if isa(v, JuMP.VariableRef)
        var_is_discrete = any(occursin(x, JuMP.name(v)) for x in ["_x", "_y", "_z"])
    else
        var_is_discrete = false
    end

    # Fix binary variables to desired values.
    _fix_indicator.(v_one, 1.0)
    _fix_indicator.(v_zero, 0.0)
    vars_relaxed = Vector{JuMP.VariableRef}([])

    if ismultinetwork(wm)
        nw_id = bound_problem.variable_to_tighten.network_index
        nw_ids_relaxed = setdiff(Set(nw_ids(wm)), Set(nw_id))
        vars_relaxed = vcat(get_all_binary_vars_at_nw!.(Ref(wm), nw_ids_relaxed)...)
        JuMP.unset_binary.(vars_relaxed)
    end

    # Optimize the variable (or affine expression) being tightened.
    JuMP.@objective(wm.model, bound_problem.sense, v)
    JuMP.optimize!(wm.model)
    termination_status = JuMP.termination_status(wm.model)

    if !(termination_status in [JuMP.INFEASIBLE, JuMP.INFEASIBLE_OR_UNBOUNDED])
        try
            # Store the candidate bound calculated from the optimization.
            candidate = JuMP.objective_bound(wm.model)
        catch
            # Otherwise, store the original bound.
            candidate = bound_problem.bound
        end
    else
        bound_problem.infeasible = true
        candidate = bound_problem.bound
    end

    if ismultinetwork(wm)
        JuMP.set_binary.(vars_relaxed)
    end

    # Update the candidate if it's for a discrete variable.
    if var_is_discrete && bound_problem.sense === JuMP.MIN_SENSE
        candidate = candidate > 0.01 ? 1.0 : candidate
    elseif var_is_discrete && bound_problem.sense === JuMP.MAX_SENSE
        candidate = candidate < 0.99 ? 0.0 : candidate
    end

    # Unfix binary variables that were fixed above.
    _unfix_indicator.(v_one), _unfix_indicator.(v_zero)

    # Return the better of the old and new bounds.
    if bound_problem.sense === JuMP.MIN_SENSE
        return max(bound_problem.bound, candidate)
    elseif bound_problem.sense === JuMP.MAX_SENSE
        return min(bound_problem.bound, candidate)
    end
end


function _set_new_bound!(bound_problem::BoundProblem, candidate::Float64)
    prec = bound_problem.precision
    num_digits = Int(ceil(-log10(prec)))
 
    if bound_problem.sense === JuMP.MIN_SENSE
        scaled_candidate = floor(candidate; digits = num_digits)

        if scaled_candidate > bound_problem.bound
            setfield!(bound_problem, :bound, scaled_candidate)
            setfield!(bound_problem, :changed, true)
        else
            setfield!(bound_problem, :changed, false)
        end
    else
        scaled_candidate = ceil(candidate; digits = num_digits)

        if scaled_candidate < bound_problem.bound
            setfield!(bound_problem, :bound, scaled_candidate)
            setfield!(bound_problem, :changed, true)
        else
            setfield!(bound_problem, :changed, false)
        end
    end
end


function _update_data_bounds!(data::Dict{String,<:Any}, bound_problems::Array{BoundProblem,1})
    data_wm = get_wm_data(data)

    for bound_problem in bound_problems
        vid = bound_problem.variable_to_tighten
        nw = string(vid.network_index)
        data_nw = ismultinetwork(data_wm) ? data_wm["nw"][nw] : data_wm
        comp = data_nw[string(vid.component_type)][string(vid.component_index)]
        comp[bound_problem.bound_name] = bound_problem.bound
    end

    # Perform data correction routines.
    correct_pipes!(data_wm)
    correct_des_pipes!(data_wm)
    correct_pumps!(data_wm)
    correct_regulators!(data_wm)
    correct_short_pipes!(data_wm)
    correct_ne_short_pipes!(data_wm)
    correct_valves!(data_wm)
    correct_nodes!(data_wm)
end


function _log_node_bound_width(data::Dict{String,<:Any}, component_type::String)
    if ismultinetwork(data)
        total_width, count = 0.0, 0

        for nw_data in values(data["nw"])
            for component in values(nw_data[component_type])
                total_width += component["head_max"] - component["head_min"]
                count += 1
            end
        end

        if count > 0
            return "Head range: $(round(total_width * inv(count); digits=3))"
        else
            return ""
        end
    elseif length(data[component_type]) > 0
        components = data[component_type]
        total_width = sum([x["head_max"] - x["head_min"] for (i, x) in components])
        return "Head range: $(round(total_width * inv(length(components)); digits=3))"
    else
        return ""
    end
end


function _log_flow_bound_width(data::Dict{String,<:Any}, component_type::String)
    if ismultinetwork(data)
        total_width, count = 0.0, 0

        for nw_data in values(data["nw"])
            for component in values(nw_data[component_type])
                total_width += component["flow_max"] - component["flow_min"]
                count += 1
            end
        end

        if count > 0
            return ", $(component_type) flow range: $(round(total_width * inv(count); digits=3))"
        else
            return ""
        end
    elseif length(data[component_type]) > 0
        components = data[component_type]
        total_width = sum([x["flow_max"] - x["flow_min"] for (i, x) in components])
        return ", $(component_type) flow range: $(round(total_width *
            inv(length(components)); digits=3))"
    else
        return ""
    end

    total_width = sum([x["flow_max"] - x["flow_min"] for (i, x) in pipes])
    return ", $(component_type) flow range:
        $(round(total_width * inv(length(pipes)); digits=3))"
end


function _log_forward_flow_bound_width(data::Dict{String,<:Any}, component_type::String)
    if ismultinetwork(data)
        total_width, count = 0.0, 0

        for nw_data in values(data["nw"])
            for component in values(nw_data[component_type])
                total_width += component["flow_max"] - component["flow_min_forward"]
                count += 1
            end
        end

        if count > 0
            return ", $(component_type) flow range: $(round(total_width * inv(count); digits=3))"
        else
            return ""
        end
    elseif length(data[component_type]) > 0
        components = data[component_type]
        total_width = sum([x["flow_max"] - x["flow_min_forward"] for (i, x) in components])
        return ", $(component_type) flow range: $(round(total_width *
            inv(length(components)); digits=3))"
    else
        return ""
    end

    total_width = sum([x["flow_max"] - x["flow_min"] for (i, x) in pipes])
    return ", $(component_type) flow range:
        $(round(total_width * inv(length(pipes)); digits=3))"
end


function _log_bound_widths(data::Dict{String,<:Any})
    wm_data = get_wm_data(data)

    message = ""
    message *= _log_node_bound_width(wm_data, "node")    
    message *= _log_flow_bound_width(wm_data, "pipe")
    message *= _log_flow_bound_width(wm_data, "des_pipe")
    message *= _log_flow_bound_width(wm_data, "short_pipe")
    message *= _log_flow_bound_width(wm_data, "ne_short_pipe")
    message *= _log_flow_bound_width(wm_data, "valve")
    message *= _log_forward_flow_bound_width(wm_data, "pump")
    message *= _log_forward_flow_bound_width(wm_data, "regulator")

    return message
end


"Remove bound problems where a discrete variable has been successfully fixed."
function _clean_bound_problems(problems::Vector{BoundProblem}, vals::Vector{Float64})
    # Initialize vectors for storing fixed binary variables.
    fixed_one_vars = Vector{_VariableIndex}([])
    fixed_zero_vars = Vector{_VariableIndex}([])

    # Enumerate over all the problems in the vector.
    for (i, problem) in enumerate(problems)
        # Get the index of the variable that was tightened.
        vid = problem.variable_to_tighten

        # If the problem fixes other variables, too, don't remove it.
        has_fixed_ones = length(problem.variables_fix_one) > 0
        has_fixed_zeros = length(problem.variables_fix_zero) > 0
        has_fixed_ones || has_fixed_zeros && (continue)

        # If the variable of interest is a model binary variable.
        if any(occursin.(["x_", "y_", "z_"], string(vid.variable_symbol)))
            # If the variable can now be fixed to one, append it to the list.
            is_fixed_to_one = vals[i] == 1.0 && problem.sense === JuMP.MIN_SENSE
            is_fixed_to_one && (append!(fixed_one_vars, [vid]))

            # If the variable can now be fixed to zero, append it to the list.
            is_fixed_to_zero = vals[i] == 0.0 && problem.sense === JuMP.MAX_SENSE
            is_fixed_to_zero && (append!(fixed_zero_vars, [vid]))
        end
    end

    # Initialize the new vector of problems, which may be reduced.
    problems_new = Vector{BoundProblem}([])
    
    for problem in problems
        # Check if the problem fixes one of the variables that has now been fixed.
        contains_fixed_one = any([x in problem.variables_fix_one for x in fixed_zero_vars])
        contains_fixed_zero = any([x in problem.variables_fix_zero for x in fixed_one_vars])

        # If the problem fixes variables that we already proved can be fixed, remove it.
        if !(contains_fixed_one || contains_fixed_zero) && !problem.infeasible
            push!(problems_new, problem)
        end
    end

    # Return the new vector of bound problems.
    return problems_new
end


"Solve a sequence of OBBT problems based on multinetwork data."
function solve_obbt_seq(data::Dict{String,<:Any}, build_method::Function, optimizer; kwargs...)
    num_time_steps = data["time_series"]["num_steps"]
    network_mn = make_multinetwork(data)

    Memento.info(_LOGGER, "[OBBT] Beginning multistep OBBT routine.")
    time_elapsed = 0.0

    for nw in 1:num_time_steps-1
        message = "[OBBT] Starting OBBT for network with multinetwork index $(nw)."
        Memento.info(_LOGGER, message)

        data_nw = deepcopy(data)
        _IM.load_timepoint!(data_nw, nw)

        if nw != 1
            map(x -> x["dispatchable"] = true, values(data_nw["tank"]))
            time_elapsed += @elapsed solve_obbt!(data_nw, build_method, optimizer; kwargs...)
            map(x -> x["dispatchable"] = false, values(data_nw["tank"]))
        else
            time_elapsed += @elapsed solve_obbt!(data_nw, build_method, optimizer; kwargs...)
        end

        network_mn["nw"][string(nw)] = data_nw

        for key in _wm_global_keys
            delete!(network_mn["nw"][string(nw)], key)
        end
    end

    time_elapsed_rounded = round(time_elapsed; digits = 2)
    Memento.info(_LOGGER, "Multistep OBBT routine completed in $(time_elapsed_rounded) seconds.")

    return network_mn
end


function solve_obbt!(
    data::Dict{String,<:Any}, build_method::Function, optimizer;
    model_type::Type = PWLRDWaterModel, time_limit::Float64 = 3600.0,
    upper_bound::Float64 = Inf, max_iter::Int = 100, relax_integrality::Bool = false,
    flow_partition_func::Function = x -> set_flow_partitions_si!(x, 1.0, 1.0e-4),
    limit_problems::Bool = false, kwargs...)
    # Print a message with relevant algorithm limit information.
    message = "[OBBT] Time limit set to $(time_limit) seconds."
    Memento.info(_LOGGER, message)

    # Execute the flow partitioning function.
    flow_partition_func(data)

    # Instantiate the bound tightening model and relax integrality, if requested.
    wms = Vector{AbstractWaterModel}(undef, Threads.nthreads())

    # Update WaterModels objects in parallel.
    Threads.@threads for i in 1:Threads.nthreads()
        wms[i] = instantiate_model(deepcopy(data), model_type, build_method)

        if upper_bound < Inf
            # Add a constraint on the objective upper bound.
            objective_expr = JuMP.objective_function(wms[i].model)
            JuMP.@constraint(wms[i].model, objective_expr <= upper_bound)
        end

        if relax_integrality
            # Relax the binary variables if requested.
            relax_all_binary_variables!(wms[i])
        end

        # Set the optimizer for the bound tightening model.
        JuMP.set_optimizer(wms[i].model, optimizer)
    end

    # Collect all bound tightening problems and update bounds in `data`.
    bound_problems = _get_bound_problems(wms[1]; limit = limit_problems)
    _update_data_bounds!(data, bound_problems) # Populate data with bounds.

    # Log mean ranges between important lower and upper bounds.
    bound_width_msg = _log_bound_widths(data)
    Memento.info(_LOGGER, "[OBBT] Initial bound widths: $(bound_width_msg).")
    Memento.info(_LOGGER, "[OBBT] Solving $(length(bound_problems)) subproblems.")

    # Instantiate termination and time logging variables.
    current_iteration = 1 # Tracks the iteration counter of the algorithm.
    time_elapsed = 0.0 # Tracks the amount of algorithm time elapsed.
    parallel_time_elapsed = 0.0 # Tracks the ideal parallel time elapsed.
    
    # Instantiate the variable that determines when OBBT is terminated.
    terminate = current_iteration >= max_iter

    while any([x.changed for x in bound_problems]) && !terminate
        # Initialize vector of bound candidates.
        vals = zeros(length(bound_problems))

        # Initialize vector of time elapsed per problem.
        parallel_times_elapsed = zeros(length(bound_problems))

        # Parallelize all bound tightening problems over available threads.
        time_elapsed += @elapsed Threads.@threads for i in 1:length(bound_problems)
            # Solve the bound tightening problem and store the candidate.
            parallel_times_elapsed = @elapsed vals[i] = _solve_bound_problem!(
                wms[Threads.threadid()], bound_problems[i])

            # Update the bound stored in the bound problem definition.
            _set_new_bound!(bound_problems[i], vals[i])
        end

        # Update the ideal parallel time elapsed from the current iteration.
        parallel_time_elapsed += maximum(parallel_times_elapsed)

        # Update all bounds within `data` using the new bounds in `bound_problems`.
        _update_data_bounds!(data, bound_problems)

        # Re-execute the flow partitioning function.
        flow_partition_func(data)

        if !terminate
            # Remove bound problems where discrete variables were fixed.
            bound_problems = _clean_bound_problems(bound_problems, vals)
        end

        # Log mean ranges between important lower and upper bounds.
        bound_width_msg = _log_bound_widths(data)
        message = "[OBBT] Iteration $(current_iteration): bound widths: $(bound_width_msg)."
        Memento.info(_LOGGER, message)
        Memento.info(_LOGGER, "[OBBT] Next iteration, Solving $(length(bound_problems)) subproblems.")

        # Check if the time limit has been exceeded and terminate, if necessary.
        time_elapsed > time_limit && ((terminate = true) && break)

        # Update the current iteration of the algorithm.
        current_iteration += 1

        # Set the termination variable if max iterations is exceeded.
        current_iteration > max_iter && (terminate = true)

        # Update WaterModels objects in parallel.
        Threads.@threads for i in 1:Threads.nthreads()
            wms[i] = instantiate_model(deepcopy(data), model_type, build_method)

            if upper_bound < Inf
                # Add a constraint on the objective upper bound.
                objective_expr = JuMP.objective_function(wms[i].model)
                JuMP.@constraint(wms[i].model, objective_expr <= upper_bound)
            end

            if relax_integrality
                # Relax the binary variables if requested.
                relax_all_binary_variables!(wms[i])
            end

            # Set the optimizer for the bound tightening model.
            JuMP.set_optimizer(wms[i].model, optimizer)
        end
    end

    time_elapsed_rounded = round(time_elapsed; digits = 2)
    parallel_time_elapsed_rounded = round(parallel_time_elapsed; digits = 2)
    Memento.info(_LOGGER, "[OBBT] Completed in $(time_elapsed_rounded) " *
        "seconds (ideal parallel time: $(parallel_time_elapsed_rounded) seconds).")
end
