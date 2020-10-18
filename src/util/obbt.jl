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
* `upper_bound`: can be used to specify a local feasible solution objective for the owf problem.
* `upper_bound_constraint`: boolean option that can be used to add an additional
   constraint to reduce the search space of each of the bound tightening
   solves. This cannot be set to `true` without specifying an upper bound.
* `use_reduced_network`: boolean option that specifies whether or not to use a reduced,
   snapshot, dispatchable version of the origin multinetwork problem for bound tightening.
"""
function run_obbt_owf!(file::String, optimizer; kwargs...)
    data = WaterModels.parse_file(file)
    return run_obbt_owf!(data, optimizer; kwargs...)
end


function _relax_model!(wm::AbstractWaterModel)
    # Further relax all binary variables in the model.
    bin_vars = filter(v -> JuMP.is_binary(v), JuMP.all_variables(wm.model))
    JuMP.unset_binary.(bin_vars) # Make all binary variables free and continuous.
    JuMP.set_lower_bound.(bin_vars, 0.0) # Lower-bound the relaxed binary variables.
    JuMP.set_upper_bound.(bin_vars, 1.0) # Upper-bound the relaxed binary variables.
end


"Translate a multinetwork dataset to a snapshot dataset with dispatchable components."
function _make_reduced_data!(ts_data::Dict{String,<:Any})
    for comp_type in keys(ts_data["time_series"])
        # If not a component type (Dict), skip parsing.
        !isa(ts_data["time_series"][comp_type], Dict) && continue

        for (i, comp) in ts_data["time_series"][comp_type]
            if comp_type == "demand"
                ts_data["demand"][i]["dispatchable"] = true
                ts_data["demand"][i]["demand_min"] = minimum(comp["flow_rate"])
                ts_data["demand"][i]["demand_max"] = maximum(comp["flow_rate"])
            elseif comp_type == "node"
                ts_data["node"][i]["h_min"] = minimum(comp["head"])
                ts_data["node"][i]["h_max"] = maximum(comp["head"])
            end
        end
    end

    for (i, tank) in ts_data["tank"]
        tank["dispatchable"] = true
    end

    for (i, reservoir) in ts_data["reservoir"]
        if "node" in keys(ts_data["time_series"])
            if string(reservoir["node"]) in keys(ts_data["time_series"]["node"])
                reservoir["dispatchable"] = true
            end
        end
    end
end


function _revert_reduced_data!(ts_data::Dict{String,<:Any})
    for comp_type in keys(ts_data["time_series"])
        # If not a component type (Dict), skip parsing.
        !isa(ts_data["time_series"][comp_type], Dict) && continue

        for (i, comp) in ts_data["time_series"][comp_type]
            if comp_type == "demand"
                ts_data["demand"][i]["dispatchable"] = false
                delete!(ts_data["demand"][i], "demand_min")
                delete!(ts_data["demand"][i], "demand_max")
            end
        end
    end

    for (i, tank) in ts_data["tank"]
        tank["dispatchable"] = false
    end

    for (i, reservoir) in ts_data["reservoir"]
        node_id = string(reservoir["node"])
        reservoir["dispatchable"] = false
        delete!(ts_data["node"][node_id], "h_min")
        delete!(ts_data["node"][node_id], "h_max")
    end
end


struct VariableIndex
    network_index::Int
    component_type::Symbol
    variable_symbol::Symbol
    component_index::Int
end


mutable struct BoundProblem
    sense::_MOI.OptimizationSense
    variable_to_tighten::VariableIndex
    variables_one::Array{VariableIndex} # Fix to one.
    variables_zero::Array{VariableIndex} # Fix to zero.
    bound_name::String # e.g., q_min, q_max, q_min_forward
    bound::Float64
    precision::Float64
    changed::Bool
end


function _get_variable_from_index(wm::AbstractWaterModel, index::VariableIndex)
    v = var(wm, index.network_index, index.variable_symbol, index.component_index)
end


function _lower_bound(wm::AbstractWaterModel, index::VariableIndex)
    v = _get_variable_from_index(wm, index)

    if typeof(v) === JuMP.VariableRef
        return JuMP.lower_bound(v)
    elseif typeof(v) == JuMP.AffExpr
        neg = collect(keys(filter(x -> x[2] < 0.0, v.terms)))
        pos = collect(keys(filter(x -> x[2] >= 0.0, v.terms)))
        return sum(JuMP.lower_bound.(pos)) - sum(JuMP.upper_bound.(neg))
    end
end


function _upper_bound(wm::AbstractWaterModel, index::VariableIndex)
    v = _get_variable_from_index(wm, index)

    if typeof(v) === JuMP.VariableRef
        return JuMP.upper_bound(v)
    elseif typeof(v) == JuMP.AffExpr
        neg = collect(keys(filter(x -> x[2] < 0.0, v.terms)))
        pos = collect(keys(filter(x -> x[2] >= 0.0, v.terms)))
        return sum(JuMP.upper_bound.(pos)) - sum(JuMP.lower_bound.(neg))
    end
end


function _create_node_bound_problems(wm::AbstractWaterModel)
    bps = Array{BoundProblem, 1}()

    for node_id in collect(ids(wm, :node))
        h = VariableIndex(wm.cnw, :node, :h, node_id)
        h_min, h_max = _lower_bound(wm, h), _upper_bound(wm, h)
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, h, [], [], "h_min", h_min, 1.0e-2, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, h, [], [], "h_max", h_max, 1.0e-2, true)])
    end

    return bps
end


function _create_pipe_bound_problems(wm::AbstractWaterModel)
    bps = Array{BoundProblem, 1}()

    for pipe_id in collect(ids(wm, :pipe))
        q = VariableIndex(wm.cnw, :pipe, :q_pipe, pipe_id)
        y = VariableIndex(wm.cnw, :pipe, :y_pipe, pipe_id)

        q_min, q_max = _lower_bound(wm, q), _upper_bound(wm, q)
        q_min_forward = get(ref(wm, q.network_index, :pipe)[pipe_id], "q_min_forward", 0.0)
        q_max_reverse = get(ref(wm, q.network_index, :pipe)[pipe_id], "q_max_reverse", 0.0)

        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [], [], "q_min", q_min, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [], [], "q_max", q_max, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [y], [], "q_min_forward", q_min_forward, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [], [y], "q_max_reverse", q_max_reverse, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, y, [], [], "y_min", 0.0, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, y, [], [], "y_max", 1.0, 1.0e-3, true)])
    end

    return bps
end


function _create_pump_bound_problems(wm::AbstractWaterModel)
    bps = Array{BoundProblem, 1}()

    for pump_id in collect(ids(wm, :pump))
        q = VariableIndex(wm.cnw, :pump, :q_pump, pump_id)
        y = VariableIndex(wm.cnw, :pump, :y_pump, pump_id)
        z = VariableIndex(wm.cnw, :pump, :z_pump, pump_id)

        q_min_forward = get(ref(wm, q.network_index, :pump)[pump_id], "q_min_forward", _FLOW_MIN)
        q_max = _upper_bound(wm, q)

        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [z, y], [], "q_min_forward", q_min_forward, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [z, y], [], "q_max", q_max, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, z, [], [], "z_min", 0.0, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, z, [], [], "z_max", 1.0, 1.0e-3, true)])
    end

    return bps
end


function _create_regulator_bound_problems(wm::AbstractWaterModel)
    bps = Array{BoundProblem, 1}()

    for regulator_id in collect(ids(wm, :regulator))
        q = VariableIndex(wm.cnw, :regulator, :q_regulator, regulator_id)
        z = VariableIndex(wm.cnw, :regulator, :z_regulator, regulator_id)

        q_min_forward = get(ref(wm, q.network_index, :regulator)[regulator_id], "q_min_forward", _FLOW_MIN)
        q_max = _upper_bound(wm, q)

        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [z], [], "q_min_forward", q_min_forward, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [z], [], "q_max", q_max, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, z, [], [], "z_min", 0.0, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, z, [], [], "z_min", 0.0, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, y, [z], [], "y_min", 0.0, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, y, [z], [], "y_max", 1.0, 1.0e-3, true)])
    end

    return bps
end


function _create_short_pipe_bound_problems(wm::AbstractWaterModel)
    bps = Array{BoundProblem, 1}()

    for short_pipe_id in collect(ids(wm, :short_pipe))
        q = VariableIndex(wm.cnw, :short_pipe, :q_short_pipe, short_pipe_id)
        y = VariableIndex(wm.cnw, :short_pipe, :y_short_pipe, short_pipe_id)

        q_min, q_max = _lower_bound(wm, q), _upper_bound(wm, q)
        q_min_forward = get(ref(wm, q.network_index, :short_pipe)[short_pipe_id], "q_min_forward", 0.0)
        q_max_reverse = get(ref(wm, q.network_index, :short_pipe)[short_pipe_id], "q_max_reverse", 0.0)

        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [], [], "q_min", q_min, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [], [], "q_max", q_max, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [y], [], "q_min_forward", q_min_forward, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [], [y], "q_max_reverse", q_max_reverse, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, y, [], [], "y_min", 0.0, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, y, [], [], "y_max", 1.0, 1.0e-3, true)])
    end

    return bps
end


function _create_valve_bound_problems(wm::AbstractWaterModel)
    bps = Array{BoundProblem, 1}()

    for valve_id in collect(ids(wm, :valve))
        q = VariableIndex(wm.cnw, :valve, :q_valve, valve_id)
        y = VariableIndex(wm.cnw, :valve, :y_valve, valve_id)
        z = VariableIndex(wm.cnw, :valve, :z_valve, valve_id)

        q_min, q_max = _lower_bound(wm, q), _upper_bound(wm, q)
        q_min_forward = get(ref(wm, q.network_index, :valve)[valve_id], "q_min_forward", 0.0)
        q_max_reverse = get(ref(wm, q.network_index, :valve)[valve_id], "q_max_reverse", 0.0)

        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [z], [], "q_min", q_min, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [z], [], "q_max", q_max, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, q, [z, y], [], "q_min_forward", q_min_forward, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, q, [z], [y], "q_max_reverse", q_max_reverse, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, y, [z], [], "y_min", 0.0, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, y, [z], [], "y_max", 1.0, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MIN_SENSE, z, [], [], "z_min", 0.0, 1.0e-3, true)])
        append!(bps, [BoundProblem(_MOI.MAX_SENSE, z, [], [], "z_max", 1.0, 1.0e-3, true)])
    end

    return bps
end


function _solve_bound_problem!(wm::AbstractWaterModel, bound_problem::BoundProblem)
    v = _get_variable_from_index(wm, bound_problem.variable_to_tighten)
    v_one = _get_variable_from_index.(Ref(wm), bound_problem.variables_one)
    v_zero = _get_variable_from_index.(Ref(wm), bound_problem.variables_zero)

    # Fix binary variables to desired values.
    JuMP.fix.(v_one, 1.0), JuMP.fix.(v_zero, 0.0)

    # Optimize the variable (or affine expression) being tightened.
    JuMP.@objective(wm.model, bound_problem.sense, v)
    JuMP.optimize!(wm.model)
    termination_status = JuMP.termination_status(wm.model)

    # Store the candidate bound calculated from the optimization.
    if termination_status in [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL]
        candidate = JuMP.objective_value(wm.model)
    end

    # Unfix binary variables that were fixed above.
    JuMP.unfix.(v_one), JuMP.unfix.(v_zero)

    # Return an optimized bound or the initial bound that was started with.
    if termination_status in [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL]
        # Get the objective value and return the better of the old and new bounds.
        if bound_problem.sense === _MOI.MIN_SENSE
            return max(bound_problem.bound, candidate)
        elseif bound_problem.sense === _MOI.MAX_SENSE
            return min(bound_problem.bound, candidate)
        end
    else
        message = "[OBBT] Optimization of $(bound_problem.variable_to_tighten) errored. Adjust tolerances."
        termination_status !== _MOI.TIME_LIMIT && Memento.warn(_LOGGER, message)
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
    total_width = sum([x["h_max"] - x["h_min"] for (i, x) in nodes])
    return "Head range: $(round(total_width * inv(length(nodes)); digits=3))"
end


function _log_pipe_bound_width(pipes::Dict{String,<:Any})
    total_width = sum([x["q_max"] - x["q_min"] for (i, x) in pipes])
    return "Pipe flow range: $(round(total_width * inv(length(pipes)); digits=3))"
end


function _log_pump_bound_width(pumps::Dict{String,<:Any})
    total_width = sum([x["q_max"] - x["q_min_forward"] for (i, x) in pumps])
    return "Pump flow range: $(round(total_width * inv(length(pumps)); digits=3))"
end


function _log_short_pipe_bound_width(short_pipes::Dict{String,<:Any})
    total_width = sum([x["q_max"] - x["q_min"] for (i, x) in short_pipes])
    return "Short pipe flow range: $(round(total_width * inv(length(short_pipes)); digits=3))"
end


function _log_regulator_bound_width(regulators::Dict{String,<:Any})
    total_width = sum([x["q_max"] - x["q_min"] for (i, x) in regulators])
    return "Regulator flow range: $(round(total_width * inv(length(regulators)); digits=3))"
end


function _log_valve_bound_width(valves::Dict{String,<:Any})
    total_width = sum([x["q_max"] - x["q_min"] for (i, x) in valves])
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


function run_obbt_owf!(data::Dict{String,<:Any}, optimizer; use_reduced_network::Bool = true,
    model_type::Type = CQRDWaterModel, time_limit::Float64 = 3600.0, upper_bound::Float64 =
    Inf, upper_bound_constraint::Bool = false, max_iter::Int = 100, improvement_tol::Float64
    = 1.0e-6, relaxed::Bool = true, precision = 1.0e-3, min_width::Float64 = 1.0e-2,
    ext::Dict{Symbol,<:Any} = Dict{Symbol,Any}(:pump_breakpoints=>3), kwargs...)
    # Print a message with relevant algorithm limit information.
    Memento.info(_LOGGER, "[OBBT] Maximum time limit for OBBT set to default value of $(time_limit) seconds.")

    if !_IM.ismultinetwork(data) && use_reduced_network
        _make_reduced_data!(data)
        build_type = WaterModels.build_wf
    elseif _IM.ismultinetwork(data) && !use_reduced_network
        build_type = WaterModels.build_mn_owf
    else
        build_type = WaterModels.build_mn_owf
        Memento.error(_LOGGER, "[OBBT] Ensure input network is the correct type.")
    end

    # Check for keyword argument inconsistencies.
    _check_obbt_options(upper_bound, upper_bound_constraint)

    # Instantiate the bound tightening model and relax integrality, if specified.
    wm = instantiate_model(data, model_type, build_type; ext=ext)
    upper_bound_constraint && _constraint_obj_bound(wm, upper_bound)
    relaxed && _relax_model!(wm) # Relax integrality, if required.

    # Set the optimizer for the bound tightening model.
    JuMP.set_optimizer(wm.model, optimizer)

    # Collect all problems.
    bound_problems = _create_bound_problems(wm)
    _update_data_bounds!(data, bound_problems) # Populate data with bounds.

    # Set up algorithm metadata.
    current_iteration = 1

    # Log widths.
    bound_width_msg = _log_bound_widths(data)
    Memento.info(_LOGGER, "[OBBT] Iteration $(current_iteration) bound widths: $(bound_width_msg).")

    while any([x.changed for x in bound_problems])
        # Obtain new candidate bounds, update bounds, and update the data.
        vals = _solve_bound_problem!.(Ref(wm), bound_problems)
        _set_new_bound!.(bound_problems, vals)
        _update_data_bounds!(data, bound_problems)

        # Set up the next optimization problem using the new bounds.
        wm = instantiate_model(data, model_type, build_type; ext=ext)
        upper_bound_constraint && _constraint_obj_bound(wm, upper_bound)
        relaxed && _relax_model!(wm)
        JuMP.set_optimizer(wm.model, optimizer)

        # Update algorithm metadata.
        current_iteration += 1

        # Log widths.
        bound_width_msg = _log_bound_widths(data)
        Memento.info(_LOGGER, "[OBBT] Iteration $(current_iteration) bound widths: $(bound_width_msg).")

        # Set the termination variable if max iterations is exceeded.
        current_iteration >= max_iter && (terminate = true)
    end

    if use_reduced_network
        _revert_reduced_data!(data)
    end
end


function _check_obbt_options(ub::Float64, ub_constraint::Bool)
    if ub_constraint && isinf(ub)
        Memento.error(_LOGGER, "[OBBT] The option \"upper_bound_constraint\" cannot be set to true without specifying an upper bound.")
    end
end
