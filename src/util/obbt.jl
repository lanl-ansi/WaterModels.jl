### Optimization-based Bound Tightening for Optimal Water Flow Problems ###
### This is built upon https://github.com/lanl-ansi/PowerModels.jl/blob/master/src/util/obbt.jl


"""
Iteratively tighten bounds on water model variables.
By default, the function uses the PWLRD relaxation for performing bound tightening.

# Example

The function can be invoked as follows:

```
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer)
solve_obbt_owf!("examples/data/epanet/van_zyl.inp", gurobi)
```

# Keyword Arguments
* `model_type`: relaxation to use for performing bound tightening.
* `max_iter`: maximum number of bound tightening iterations to perform.
* `time_limit`: maximum amount of time (seconds) for the bound tightening algorithm.
* `upper_bound`: can be used to specify a local feasible solution objective for the problem.
* `upper_bound_constraint`: boolean option that can be used to add an additional
   constraint to reduce the search space of each of the bound tightening
   solves. This cannot be set to `true` without specifying an upper bound.
* `use_relaxed_network`: boolean option that specifies whether or not to use a relaxed,
   snapshot, dispatchable version of the origin multinetwork problem for bound tightening.
"""
function solve_obbt_owf!(file::String, optimizer; kwargs...)
    data = WaterModels.parse_file(file)
    return solve_obbt_owf!(data, optimizer; kwargs...)
end


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
    _fix_indicator.(v_one, 1.0), _fix_indicator.(v_zero, 0.0)
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

    if ismultinetwork(wm)
        JuMP.set_binary.(vars_relaxed)
    end

    # Store the candidate bound calculated from the optimization.
    if termination_status in [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL]
        candidate = JuMP.objective_value(wm.model)

        if var_is_discrete && bound_problem.sense === _MOI.MIN_SENSE
            candidate = candidate > 0.01 ? 1.0 : candidate
        elseif var_is_discrete && bound_problem.sense === _MOI.MAX_SENSE
            candidate = candidate < 0.99 ? 0.0 : candidate
        end
    end

    # Unfix binary variables that were fixed above.
    _unfix_indicator.(v_one), _unfix_indicator.(v_zero)

    # Return an optimized bound or the initial bound that was started with.
    if termination_status in [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL]
        # Get the objective value and return the better of the old and new bounds.
        if bound_problem.sense === _MOI.MIN_SENSE
            return max(bound_problem.bound, candidate)
        elseif bound_problem.sense === _MOI.MAX_SENSE
            return min(bound_problem.bound, candidate)
        end
    else
        return bound_problem.bound # Optimization was not successful. Return the starting bound.
    end
end


function _set_new_bound!(bound_problem::BoundProblem, candidate::Float64)
    prec = bound_problem.precision

    if bound_problem.sense === _MOI.MIN_SENSE
        scaled_candidate = floor(inv(prec) * candidate) * prec

        if scaled_candidate > bound_problem.bound
            setfield!(bound_problem, :bound, scaled_candidate)
            setfield!(bound_problem, :changed, true)
        else
            setfield!(bound_problem, :changed, false)
        end
    else
        scaled_candidate = ceil(inv(prec) * candidate) * prec

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
    message *= _log_flow_bound_width(wm_data, "short_pipe")
    message *= _log_flow_bound_width(wm_data, "valve")
    message *= _log_forward_flow_bound_width(wm_data, "pump")
    message *= _log_forward_flow_bound_width(wm_data, "regulator")

    return message
end


"Remove bound problems where a discrete variable has been successfully fixed."
function _clean_bound_problems!(problems::Vector{BoundProblem}, vals::Vector{Float64})
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
            is_fixed_to_one = vals[i] == 1.0 && problem.sense === _MOI.MIN_SENSE
            is_fixed_to_one && (append!(fixed_one_vars, [vid]))

            # If the variable can now be fixed to zero, append it to the list.
            is_fixed_to_zero = vals[i] == 0.0 && problem.sense === _MOI.MAX_SENSE
            is_fixed_to_zero && (append!(fixed_zero_vars, [vid]))
        end
    end

    for problem in problems
        # Check if the problem fixes one of the variables that has now been fixed.
        contains_fixed_one = any([x in problem.variables_fix_one for x in fixed_zero_vars])
        contains_fixed_zero = any([x in problem.variables_fix_zero for x in fixed_one_vars])

        # If the problem fixes variables that we already proved can be fixed, remove it.
        if contains_fixed_one || contains_fixed_zero
            problems = Vector{BoundProblem}(setdiff!(problems, [problem]))
        end
    end
end


function solve_obbt_owf!(
    data::Dict{String,<:Any}, optimizer; use_relaxed_network::Bool = true,
    model_type::Type = PWLRDWaterModel, time_limit::Float64 = 3600.0,
    upper_bound::Float64 = Inf, upper_bound_constraint::Bool = false,
    max_iter::Int = 100, solve_relaxed::Bool = true,
    flow_partition_func::Function = x -> set_flow_partitions_si!(x, 1.0, 1.0e-4),
    limit_problems::Bool = false, kwargs...)
    # Print a message with relevant algorithm limit information.
    message = "[OBBT] Maximum time limit set to default value of $(time_limit) seconds."
    Memento.info(_LOGGER, message)

    # Relax the network (e.g., make nodal components dispatchable) if requested.
    use_relaxed_network && relax_network!(data)

    # Set the problem specification that will be used for bound tightening.
    build_type = _IM.ismultinetwork(get_wm_data(data)) ? build_mn_owf : build_wf

    # Check for keyword argument inconsistencies.
    _check_obbt_options(upper_bound, upper_bound_constraint)

    # Instantiate the bound tightening model and relax integrality, if requested.
    wms = [instantiate_model(data, model_type, build_type) for i in 1:Threads.nthreads()]

    if upper_bound_constraint
        map(x -> _constraint_obj_bound(x, upper_bound), wms)
    end

    if solve_relaxed
        # Relax the binary variables if requested.
        map(x -> relax_all_binary_variables!(x), wms)
    end

    # Set the optimizer for the bound tightening model.
    map(x -> JuMP.set_optimizer(x.model, optimizer), wms)

    # Collect all problems.
    bound_problems = _get_bound_problems(wms[1]; limit = limit_problems)
    _update_data_bounds!(data, bound_problems) # Populate data with bounds.

    # Log widths.
    bound_width_msg = _log_bound_widths(data)
    Memento.info(_LOGGER, "[OBBT] Initial bound widths: $(bound_width_msg).")
    terminate, time_elapsed = false, 0.0

    # Set up algorithm metadata.
    current_iteration = 1
    terminate = current_iteration >= max_iter

    while any([x.changed for x in bound_problems]) && !terminate
        # Obtain new candidate bounds, update bounds, and update the data.
        vals = zeros(length(bound_problems))

        time_elapsed += @elapsed Threads.@threads for i in 1:length(bound_problems)
            vals[i] = _solve_bound_problem!(wms[Threads.threadid()], bound_problems[i])
            _set_new_bound!(bound_problems[i], vals[i])
        end

        time_elapsed > time_limit && ((terminate = true) && break)
        _update_data_bounds!(data, bound_problems)
        flow_partition_func(data)
        !terminate && _clean_bound_problems!(bound_problems, vals)

        # Log widths.
        bound_width_msg = _log_bound_widths(data)
        message = "[OBBT] Iteration $(current_iteration) bound widths: $(bound_width_msg)."
        Memento.info(_LOGGER, message)

        # Update algorithm metadata.
        current_iteration += 1

        # Set up the next optimization problem using the new bounds.
        wms = [instantiate_model(data, model_type, build_type) for i in 1:Threads.nthreads()]

        if upper_bound_constraint
            map(x -> _constraint_obj_bound(x, upper_bound), wms)
        end

        if solve_relaxed
            # Relax the binary variables if requested.
            map(x -> relax_all_binary_variables!(x), wms)
        end

        # Set the optimizer for the bound tightening model.
        map(x -> JuMP.set_optimizer(x.model, optimizer), wms)

        # Set the termination variable if max iterations is exceeded.
        current_iteration >= max_iter && (terminate = true)
    end

    if use_relaxed_network
        _fix_demands!(data)
        _fix_tanks!(data)
        _fix_reservoirs!(data)
    end

    time_elapsed_rounded = round(time_elapsed; digits = 2)
    Memento.info(_LOGGER, "[OBBT] Completed in $(time_elapsed_rounded) seconds.")
end


function _check_obbt_options(ub::Float64, ub_constraint::Bool)
    if ub_constraint && isinf(ub)
        message = "[OBBT] The option \"upper_bound_constraint\" cannot " *
            "be set to true without specifying an upper bound."
        Memento.error(_LOGGER, message)
    end
end
