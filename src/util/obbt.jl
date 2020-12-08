### Optimization-based Bound Tightening for Optimal Water Flow Problems ###
### This is built upon https://github.com/lanl-ansi/PowerModels.jl/blob/master/src/util/obbt.jl


"""
Iteratively tighten bounds on head and flow variables.

The function can be invoked on any convex relaxation which explicitly has these variables.
By default, the function uses the CRD relaxation for performing bound tightening.

# Example

The function can be invoked as follows:

```
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer)
run_obbt_owf!("examples/data/epanet/van_zyl.inp", ipopt)
```

# Keyword Arguments
* `model_type`: relaxation to use for performing bound tightening.
* `max_iter`: maximum number of bound tightening iterations to perform.
* `min_width`: domain width beyond which bound tightening is not performed.
* `precision`: decimal precision to round the tightened bounds to.
* `time_limit`: maximum amount of time (seconds) for the bound tightening algorithm.
* `upper_bound`: can be used to specify a local feasible solution objective for the problem.
* `upper_bound_constraint`: boolean option that can be used to add an additional
   constraint to reduce the search space of each of the bound tightening
   solves. This cannot be set to `true` without specifying an upper bound.
* `use_relaxed_network`: boolean option that specifies whether or not to use a relaxed,
   snapshot, dispatchable version of the origin multinetwork problem for bound tightening.
"""
function run_obbt_owf!(file::String, optimizer; kwargs...)
    data = WaterModels.parse_file(file)
    return run_obbt_owf!(data, optimizer; kwargs...)
end


function _revert_relaxed_data!(ts_data::Dict{String,<:Any})
    for comp_type in keys(ts_data["time_series"])
        # If not a component type (Dict), skip parsing.
        !isa(ts_data["time_series"][comp_type], Dict) && continue

        for (i, comp) in ts_data["time_series"][comp_type]
            if comp_type == "demand"
                ts_data["demand"][i]["dispatchable"] = false
                #delete!(ts_data["demand"][i], "flow_min")
                #delete!(ts_data["demand"][i], "flow_max")
            end
        end
    end

    for (i, tank) in ts_data["tank"]
        tank["dispatchable"] = false
    end

    for (i, reservoir) in ts_data["reservoir"]
        #node_id = string(reservoir["node"])
        reservoir["dispatchable"] = false
        #delete!(ts_data["node"][node_id], "head_min")
        #delete!(ts_data["node"][node_id], "head_max")
    end
end


mutable struct BoundProblem
    sense::_MOI.OptimizationSense
    variable_to_tighten::_VariableIndex
    variables_fix_one::Array{_VariableIndex} # Fix to one.
    variables_fix_zero::Array{_VariableIndex} # Fix to zero.
    bound_name::String # e.g., flow_min, flow_max, flow_min_forward
    bound::Float64
    precision::Float64
    changed::Bool
end


function _create_node_bound_problems(wm::AbstractWaterModel)
    bps = Array{BoundProblem, 1}()

    for node_id in collect(ids(wm, :node))
        h = _VariableIndex(wm.cnw, :node, :h, node_id)
        h_min, h_max = _get_lower_bound_from_index(wm, h), _get_upper_bound_from_index(wm, h)
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, h, [], [], "head_min", h_min, 1.0e-2, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, h, [], [], "head_max", h_max, 1.0e-2, true)])
    end

    return bps
end


function _create_pipe_bound_problems(wm::AbstractWaterModel)
    bps = Array{BoundProblem, 1}()

    for pipe_id in collect(ids(wm, :pipe))
        q = _VariableIndex(wm.cnw, :pipe, :q_pipe, pipe_id)
        y = _VariableIndex(wm.cnw, :pipe, :y_pipe, pipe_id)

        flow_min, flow_max = _get_lower_bound_from_index(wm, q), _get_upper_bound_from_index(wm, q)
        flow_min_forward = get(ref(wm, q.network_index, :pipe)[pipe_id], "flow_min_forward", 0.0)
        flow_max_reverse = get(ref(wm, q.network_index, :pipe)[pipe_id], "flow_max_reverse", 0.0)

        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [], [], "flow_min", flow_min, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [], [], "flow_max", flow_max, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [y], [], "flow_min_forward", flow_min_forward, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [], [y], "flow_max_reverse", flow_max_reverse, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, y, [], [], "y_min", 0.0, 1.0e-2, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, y, [], [], "y_max", 1.0, 1.0e-2, true)])
    end

    return bps
end


function _create_pump_bound_problems(wm::AbstractWaterModel)
    bps = Array{BoundProblem, 1}()

    for pump_id in collect(ids(wm, :pump))
        q = _VariableIndex(wm.cnw, :pump, :q_pump, pump_id)
        y = _VariableIndex(wm.cnw, :pump, :y_pump, pump_id)
        z = _VariableIndex(wm.cnw, :pump, :z_pump, pump_id)

        flow_min_forward = get(ref(wm, q.network_index, :pump)[pump_id], "flow_min_forward", _FLOW_MIN)
        flow_max = _get_upper_bound_from_index(wm, q)

        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [z, y], [], "flow_min_forward", flow_min_forward, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [z, y], [], "flow_max", flow_max, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, z, [], [], "z_min", 0.0, 1.0e-2, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, z, [], [], "z_max", 1.0, 1.0e-2, true)])
    end

    return bps
end


function _create_regulator_bound_problems(wm::AbstractWaterModel)
    bps = Array{BoundProblem, 1}()

    for regulator_id in collect(ids(wm, :regulator))
        q = _VariableIndex(wm.cnw, :regulator, :q_regulator, regulator_id)
        y = _VariableIndex(wm.cnw, :regulator, :y_regulator, regulator_id)
        z = _VariableIndex(wm.cnw, :regulator, :z_regulator, regulator_id)

        flow_min_forward = get(ref(wm, q.network_index, :regulator)[regulator_id], "flow_min_forward", _FLOW_MIN)
        flow_max = _get_upper_bound_from_index(wm, q)

        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [z], [], "flow_min_forward", flow_min_forward, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [z], [], "flow_max", flow_max, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, z, [], [], "z_min", 0.0, 1.0e-2, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, z, [], [], "z_max", 1.0, 1.0e-2, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, y, [z], [], "y_min", 0.0, 1.0e-2, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, y, [z], [], "y_max", 1.0, 1.0e-2, true)])
    end

    return bps
end


function _create_short_pipe_bound_problems(wm::AbstractWaterModel)
    bps = Array{BoundProblem, 1}()

    for short_pipe_id in collect(ids(wm, :short_pipe))
        q = _VariableIndex(wm.cnw, :short_pipe, :q_short_pipe, short_pipe_id)
        y = _VariableIndex(wm.cnw, :short_pipe, :y_short_pipe, short_pipe_id)

        flow_min, flow_max = _get_lower_bound_from_index(wm, q), _get_upper_bound_from_index(wm, q)
        flow_min_forward = get(ref(wm, q.network_index, :short_pipe)[short_pipe_id], "flow_min_forward", 0.0)
        flow_max_reverse = get(ref(wm, q.network_index, :short_pipe)[short_pipe_id], "flow_max_reverse", 0.0)

        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [], [], "flow_min", flow_min, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [], [], "flow_max", flow_max, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [y], [], "flow_min_forward", flow_min_forward, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [], [y], "flow_max_reverse", flow_max_reverse, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, y, [], [], "y_min", 0.0, 1.0e-2, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, y, [], [], "y_max", 1.0, 1.0e-2, true)])
    end

    return bps
end


function _create_valve_bound_problems(wm::AbstractWaterModel)
    bps = Array{BoundProblem, 1}()

    for valve_id in collect(ids(wm, :valve))
        q = _VariableIndex(wm.cnw, :valve, :q_valve, valve_id)
        y = _VariableIndex(wm.cnw, :valve, :y_valve, valve_id)
        z = _VariableIndex(wm.cnw, :valve, :z_valve, valve_id)

        flow_min, flow_max = _get_lower_bound_from_index(wm, q), _get_upper_bound_from_index(wm, q)
        flow_min_forward = get(ref(wm, q.network_index, :valve)[valve_id], "flow_min_forward", 0.0)
        flow_max_reverse = get(ref(wm, q.network_index, :valve)[valve_id], "flow_max_reverse", 0.0)

        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [z], [], "flow_min", flow_min, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [z], [], "flow_max", flow_max, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [z, y], [], "flow_min_forward", flow_min_forward, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [z], [y], "flow_max_reverse", flow_max_reverse, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, y, [z], [], "y_min", 0.0, 1.0e-2, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, y, [z], [], "y_max", 1.0, 1.0e-2, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, z, [], [], "z_min", 0.0, 1.0e-2, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, z, [], [], "z_max", 1.0, 1.0e-2, true)])
    end

    return bps
end


function _fix_indicator(v::JuMP.VariableRef, value::Float64)
    if JuMP.is_binary(v)
        JuMP.fix(v, value)
    else
        JuMP.fix(v, value; force=true)
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

    # Fix binary variables to desired values.
    _fix_indicator.(v_one, 1.0), _fix_indicator.(v_zero, 0.0)

    # Optimize the variable (or affine expression) being tightened.
    JuMP.@objective(wm.model, bound_problem.sense, v)
    JuMP.optimize!(wm.model)
    termination_status = JuMP.termination_status(wm.model)

    # Store the candidate bound calculated from the optimization.
    if termination_status in [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL]
        candidate = JuMP.objective_value(wm.model)
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
    for bound_problem in bound_problems
        vid = bound_problem.variable_to_tighten
        comp = data[string(vid.component_type)][string(vid.component_index)]
        comp[bound_problem.bound_name] = bound_problem.bound
    end
end


function _create_bound_problems(wm::AbstractWaterModel)
    # Create the sets of bound-tightening problems.
    bps_node = _create_node_bound_problems(wm)
    bps_pipe = _create_pipe_bound_problems(wm)
    bps_pump = _create_pump_bound_problems(wm)
    bps_regulator = _create_regulator_bound_problems(wm)
    bps_short_pipe = _create_short_pipe_bound_problems(wm)
    bps_valve = _create_valve_bound_problems(wm)

    return vcat(bps_pump, bps_valve, bps_regulator, bps_pipe, bps_short_pipe, bps_node)
end


function _log_node_bound_width(nodes::Dict{String,<:Any})
    total_width = sum([x["head_max"] - x["head_min"] for (i, x) in nodes])
    return "Head range: $(round(total_width * inv(length(nodes)); digits=3))"
end


function _log_pipe_bound_width(pipes::Dict{String,<:Any})
    total_width = sum([x["flow_max"] - x["flow_min"] for (i, x) in pipes])
    return "Pipe flow range: $(round(total_width * inv(length(pipes)); digits=3))"
end


function _log_pump_bound_width(pumps::Dict{String,<:Any})
    total_width = sum([x["flow_max"] - x["flow_min_forward"] for (i, x) in pumps])
    return "Pump flow range: $(round(total_width * inv(length(pumps)); digits=3))"
end


function _log_short_pipe_bound_width(short_pipes::Dict{String,<:Any})
    total_width = sum([x["flow_max"] - x["flow_min"] for (i, x) in short_pipes])
    return "Short pipe flow range: $(round(total_width * inv(length(short_pipes)); digits=3))"
end


function _log_regulator_bound_width(regulators::Dict{String,<:Any})
    total_width = sum([x["flow_max"] - x["flow_min_forward"] for (i, x) in regulators])
    return "Regulator flow range: $(round(total_width * inv(length(regulators)); digits=3))"
end


function _log_valve_bound_width(valves::Dict{String,<:Any})
    total_width = sum([x["flow_max"] - x["flow_min"] for (i, x) in valves])
    return "Valve flow range: $(round(total_width * inv(length(valves)); digits=3))"
end


function _log_bound_widths(data::Dict{String,<:Any})
    message = ""

    message *= length(data["node"]) > 0 ? _log_node_bound_width(data["node"]) : ""
    message *= length(data["pipe"]) > 0 ? ", " * _log_pipe_bound_width(data["pipe"]) : ""
    message *= length(data["pump"]) > 0 ? ", " * _log_pump_bound_width(data["pump"]) : ""
    message *= length(data["short_pipe"]) > 0 ? ", " * _log_short_pipe_bound_width(data["short_pipe"]) : ""
    message *= length(data["regulator"]) > 0 ? ", " * _log_regulator_bound_width(data["regulator"]) : ""
    message *= length(data["valve"]) > 0 ? ", " * _log_valve_bound_width(data["valve"]) : ""

    return message
end


function _clean_bound_problems!(problems::Array{BoundProblem, 1}, vals::Array{Float64,1})
    fixed_one_vars = Array{_VariableIndex, 1}([])
    fixed_zero_vars = Array{_VariableIndex, 1}([])

    for (i, problem) in enumerate(problems)
        vid = problem.variable_to_tighten

        if length(problem.variables_fix_one) > 0 || length(problem.variables_fix_zero) > 0
            continue
        end

        if any(occursin.(["y_", "z_"], string(vid.variable_symbol)))
            is_fixed_to_one = vals[i] == 1.0 && problem.sense === _MOI.MIN_SENSE
            is_fixed_to_one && (append!(fixed_one_vars, [vid]))
            is_fixed_to_zero = vals[i] == 0.0 && problem.sense === _MOI.MAX_SENSE
            is_fixed_to_zero && (append!(fixed_zero_vars, [vid]))
        end
    end

    for (i, problem) in enumerate(problems)
        contains_fixed_one = any([x in problem.variables_fix_one for x in fixed_zero_vars])
        contains_fixed_zero = any([x in problem.variables_fix_zero for x in fixed_one_vars])

        if contains_fixed_one || contains_fixed_zero
            problems = setdiff!(problems, [problem])
        end
    end
end


function run_obbt_owf!(data::Dict{String,<:Any}, optimizer; use_relaxed_network::Bool = true,
    model_type::Type = CQRDWaterModel, time_limit::Float64 = 3600.0, upper_bound::Float64 =
    Inf, upper_bound_constraint::Bool = false, max_iter::Int = 100, improvement_tol::Float64
    = 1.0e-6, solve_relaxed::Bool = true, precision = 1.0e-2, min_width::Float64 = 1.0e-2,
    ext::Dict{Symbol,<:Any} = Dict{Symbol,Any}(:pump_breakpoints=>3), kwargs...)
    # Print a message with relevant algorithm limit information.
    Memento.info(_LOGGER, "[OBBT] Maximum time limit for OBBT set to default value of $(time_limit) seconds.")

    # Relax the network (i.e., make nodal components dispatchable) if requested.
    use_relaxed_network && _relax_network!(data)

    # Set the problem specification that will be used for bound tightening.
    build_type = _IM.ismultinetwork(data) ? build_wf : build_mn_wf

    # Check for keyword argument inconsistencies.
    _check_obbt_options(upper_bound, upper_bound_constraint)

    # Instantiate the bound tightening model and relax integrality, if requested.
    wm = instantiate_model(data, model_type, build_type; ext=ext)
    upper_bound_constraint && _constraint_obj_bound(wm, upper_bound)
    solve_relaxed && relax_all_binary_variables!(wm)

    # Set the optimizer for the bound tightening model.
    JuMP.set_optimizer(wm.model, optimizer)

    # Collect all problems.
    bound_problems = _create_bound_problems(wm)
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

        for (i, bound_problem) in enumerate(bound_problems)
            time_elapsed += @elapsed vals[i] = _solve_bound_problem!(wm, bound_problem)
            _set_new_bound!(bound_problem, vals[i])
            time_elapsed > time_limit && ((terminate = true) && break)
        end

        _update_data_bounds!(data, bound_problems)
        !terminate && _clean_bound_problems!(bound_problems, vals)

        # Log widths.
        bound_width_msg = _log_bound_widths(data)
        Memento.info(_LOGGER, "[OBBT] Iteration $(current_iteration) bound widths: $(bound_width_msg).")

        # Update algorithm metadata.
        current_iteration += 1

        # Set up the next optimization problem using the new bounds.
        wm = instantiate_model(data, model_type, build_type; ext=ext)
        upper_bound_constraint && _constraint_obj_bound(wm, upper_bound)
        solve_relaxed && relax_all_binary_variables!(wm)
        JuMP.set_optimizer(wm.model, optimizer)

        # Set the termination variable if max iterations is exceeded.
        current_iteration >= max_iter && (terminate = true)
    end

    if use_relaxed_network
        _revert_relaxed_data!(data)
    end
end


function _check_obbt_options(ub::Float64, ub_constraint::Bool)
    if ub_constraint && isinf(ub)
        Memento.error(_LOGGER, "[OBBT] The option \"upper_bound_constraint\" cannot be set to true without specifying an upper bound.")
    end
end
