### Optimization-based Bound Tightening for Optimal Water Flow Problems ###
### This is built upon https://github.com/lanl-ansi/WaterModels.jl/blob/master/src/util/obbt.jl

"""
Iteratively tighten bounds on head and flow variables.

The function can be invoked on any convex relaxation which explicitly has these variables.
By default, the function uses the MICP-R relaxation for performing bound tightening.
Interested readers are referred to the paper "Strengthening Convex Relaxations with Bound
Tightening for Water Network Optimization."

# Example

The function can be invoked as follows:

```
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer)
data, stats = run_obbt_owf!("examples/data/epanet/van_zyl.inp", ipopt)
```

`data` contains the parsed network data with tightened bounds.
`stats` contains information output from the bound tightening algorithm.

# Keyword Arguments
* `model_type`: relaxation to use for performing bound tightening.
* `max_iter`: maximum number of bound tightening iterations to perform.
* `time_limit`: maximum amount of time (seconds) for the bound tightening algorithm.
* `upper_bound`: can be used to specify a local feasible solution objective for the owf problem.
* `upper_bound_constraint`: boolean option that can be used to add an additional
   constraint to reduce the search space of each of the bound tightening
   solves. This cannot be set to `true` without specifying an upper bound.
* `rel_gap_tol`: tolerance used to terminate the algorithm when the objective
   value of the relaxation is close to the upper bound specified using the
   `upper_bound` keyword.
* `min_width`: domain beyond which bound tightening is not performed.
* `termination`: Bound-tightening algorithm terminates if the improvement in
   the average or maximum bound improvement, specified using either the
   `termination = :avg` or the `termination = :max` option, is less than
   `improvement_tol`.
* `precision`: number of decimal digits to round the tightened bounds to.
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


function _get_node_bound_dict(wm::AbstractWaterModel, nw::Int)
    return Dict(string(i) => Dict("h_min"=>JuMP.lower_bound(var(wm, nw, :h, i)),
        "h_max"=>JuMP.upper_bound(var(wm, nw, :h, i))) for i in ids(wm, nw, :node))
end


function _get_edge_bound_dict(wm::AbstractWaterModel, nw::Int, comp::Symbol)
    qp_sym, qn_sym = Symbol("qp_" * string(comp)), Symbol("qn_" * string(comp))
    qp, qn = var(wm, nw, qp_sym), var(wm, nw, qn_sym)
    return Dict(string(a) => Dict("q_min"=>min(0.0, -JuMP.upper_bound(qn[a])),
        "q_max"=>max(0.0, JuMP.upper_bound(qp[a]))) for a in ids(wm, nw, comp))
end


function _create_modifications_reduced(wm::AbstractWaterModel)
    data = Dict{String,Any}("per_unit"=>false)
    data["node"] = _get_node_bound_dict(wm, wm.cnw)
    data["pipe"] = Dict{String,Any}()

    for type in [:pipe, :shutoff_valve, :check_valve]
        qp_sym, qn_sym = Symbol("qp_" * string(type)), Symbol("qn_" * string(type))
        tmp_data = _get_edge_bound_dict(wm, wm.cnw, type)
        data["pipe"] = merge(data["pipe"], tmp_data)
    end

    for type in [:pressure_reducing_valve, :pump]
        qp_sym, qn_sym = Symbol("qp_" * string(type)), Symbol("qn_" * string(type))
        data[string(type)] = _get_edge_bound_dict(wm, wm.cnw, type)
    end

    return data
end


function _create_modifications_mn(wm::AbstractWaterModel)
    data = Dict{String,Any}("nw"=>Dict{String,Any}(), "per_unit"=>false)

    for nw in sort(collect(nw_ids(wm)))
        dnw = data["nw"][string(nw)] = Dict{String,Any}()
        dnw["node"] = _get_node_bound_dict(wm, nw)
        dnw["pipe"] = Dict{String,Any}()

        for type in [:pipe, :shutoff_valve, :check_valve]
            qp_sym, qn_sym = Symbol("qp_" * string(type)), Symbol("qn_" * string(type))
            tmp_dnw = _get_edge_bound_dict(wm, nw, type)
            dnw["pipe"] = merge(dnw["pipe"], tmp_dnw)
        end

        for type in [:pressure_reducing_valve, :pump]
            qp_sym, qn_sym = Symbol("qp_" * string(type)), Symbol("qn_" * string(type))
            dnw[string(type)] = _get_edge_bound_dict(wm, nw, type)
        end
    end

    return data
end


function _get_head_index_set(wm::AbstractWaterModel; width::Float64=1.0e-6, prec::Int=10)
    return vcat([[(nw, :node, :h, i, width, prec) for i in ids(wm, nw, :node)]
        for nw in sort(collect(nw_ids(wm)))]...)
end


function _get_flow_index_set(wm::AbstractWaterModel; width::Float64=1.0e-8, prec::Int=12)
    types = [:pipe, :shutoff_valve, :check_valve, :pressure_reducing_valve, :pump]
    return vcat([vcat([[(nw, type, Symbol("q_" * string(type)), a, width, prec) for a in
        ids(wm, nw, type)] for nw in sort(collect(nw_ids(wm)))]...) for type in types]...)
end


function _get_existing_bounds(data::Dict{String,<:Any}, index::Tuple)
    nw_str, comp_str, index_str = string(index[1]), string(index[2]), string(index[4])
    comp_str = comp_str in ["check_valve", "shutoff_valve"] ? "pipe" : comp_str
    var_str = occursin("q", string(index[3])) ? "q" : string(index[3])

    if "nw" in keys(data)
        data_index = data["nw"][nw_str][comp_str][index_str]
    else
        data_index = data[comp_str][index_str]
    end

    return data_index[var_str * "_min"], data_index[var_str * "_max"]
end


function _update_modifications!(data::Dict{String,Any}, var_index_set::Array, bounds::Array)
    # Loop over variables and tighten bounds.
    for (i, id) in enumerate(var_index_set)
        nw_str, comp_str, index_str = string(id[1]), string(id[2]), string(id[4])
        var_str = occursin("q", string(id[3])) ? "q" : string(id[3])
        comp_c = comp_str in ["shutoff_valve", "check_valve"] ? "pipe" : comp_str

        if _IM.ismultinetwork(data)
            data["nw"][nw_str][comp_c][index_str][var_str * "_min"] = bounds[i][1]
            data["nw"][nw_str][comp_c][index_str][var_str * "_max"] = bounds[i][2]
        else
            data[comp_c][index_str][var_str * "_min"] = bounds[i][1]
            data[comp_c][index_str][var_str * "_max"] = bounds[i][2]
        end
    end
end


function _get_average_widths(data::Dict{String,<:Any})
    h = sum(x["h_max"] - x["h_min"] for (i, x) in data["node"])
    message = "[OBBT] Average bound widths: h -> $(h * inv(length(data["node"]))), "
    avg_vals = [h * inv(length(data["node"]))]

    for type in ["pipe", "pressure_reducing_valve", "pump"]
        if length(data[type]) > 0
            q_length = length(data[type])
            q = sum(x["q_max"] - x["q_min"] for (i, x) in data[type])
            message *= "q_$(type) -> $(q * inv(q_length)), "
            avg_vals = vcat(avg_vals, q * inv(q_length))
        end
    end

    Memento.info(_LOGGER, message[1:end-2] * ".")
    return avg_vals
end


function _get_average_widths_mn(data::Dict{String,<:Any})
    h = sum(sum(x["h_max"] - x["h_min"] for (i, x) in nw["node"]) for (n, nw) in data["nw"])
    h_length = sum(length(nw["node"]) for (n, nw) in data["nw"])
    message = "[OBBT] Average bound widths: h -> $(h * inv(h_length)), "
    avg_vals = [h * inv(length(data["node"]))]

    for type in ["pipe", "pressure_reducing_valve", "pump"]
        if sum(length(nw[type]) for (n, nw) in data["nw"]) > 0
            q_length = sum(length(nw[type]) for (n, nw) in data["nw"])
            q = sum(sum(x["q_max"] - x["q_min"] for (i, x) in nw[type]) for (n, nw) in data["nw"])
            message *= "q_$(type) -> $(q * inv(q_length)), "
            avg_vals = vcat(avg_vals, q * inv(q_length))
        end
    end

    Memento.info(_LOGGER, message[1:end-2] * ".")
    return avg_vals
end


"Translate a multinetwork dataset to a snapshot dataset with dispatchable components."
function _make_reduced_data!(ts_data::Dict{String,<:Any})
    for comp_type in keys(ts_data["time_series"])
        # If not a component type (Dict), skip parsing.
        !isa(ts_data["time_series"][comp_type], Dict) && continue

        for (i, comp) in ts_data["time_series"][comp_type]
            if comp_type == "junction"
                ts_data["junction"][i]["dispatchable"] = true
                ts_data["junction"][i]["demand_min"] = minimum(comp["demand"])
                ts_data["junction"][i]["demand_max"] = maximum(comp["demand"])
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
        reservoir["dispatchable"] = true
    end
end


function _revert_reduced_data!(ts_data::Dict{String,<:Any})
    for comp_type in keys(ts_data["time_series"])
        # If not a component type (Dict), skip parsing.
        !isa(ts_data["time_series"][comp_type], Dict) && continue

        for (i, comp) in ts_data["time_series"][comp_type]
            if comp_type == "junction"
                ts_data["junction"][i]["dispatchable"] = false
                delete!(ts_data["junction"][i], "demand_min")
                delete!(ts_data["junction"][i], "demand_max")
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


function run_obbt_owf!(data::Dict{String,<:Any}, optimizer; use_reduced_network::Bool = true,
    model_type::Type = MICPRWaterModel, time_limit::Float64 = 3600.0, upper_bound::Float64 =
    Inf, upper_bound_constraint::Bool = false, rel_gap_tol::Float64 = Inf, max_iter::Int = 100,
    improvement_tol::Float64 = 1.0e-6, termination::Symbol = :avg, relaxed::Bool = true,
    ext::Dict{Symbol,<:Any} = Dict{Symbol,Any}(:pump_breakpoints=>5), kwargs...)
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
    _check_obbt_options(upper_bound, rel_gap_tol, upper_bound_constraint)

    # Instantiate the bound tightening model and relax integrality, if specified.
    bt = instantiate_model(data, model_type, build_type; ext=ext)
    upper_bound_constraint && _constraint_obj_bound(bt, upper_bound)
    relaxed && _relax_model!(bt) # Relax integrality, if required.

    # Build the dictionary and sets that will store to the network.
    modifications = use_reduced_network ? _create_modifications_reduced(bt) : _create_modifications_mn(bt)
    var_index_set = vcat(_get_head_index_set(bt), _get_flow_index_set(bt))
    #var_index_set = _get_flow_index_set(bt)
    bounds = [_get_existing_bounds(modifications, vid) for vid in var_index_set]

    # Get the initial average bound widths.
    if use_reduced_network
        avg_widths_initial = _get_average_widths(modifications)
    else
        avg_widths_initial = _get_average_widths_mn(modifications)
    end

    # Initialize algorithm termination metadata.
    current_iteration, time_elapsed = 1, 0.0
    parallel_time_elapsed, terminate = 0.0, false

    # Set the optimizer for the bound tightening model.
    JuMP.set_optimizer(bt.model, optimizer)

    while !terminate # Algorithmic loop.
        Memento.info(_LOGGER, "[OBBT] Starting iteration $(current_iteration).")

        # Loop over variables and tighten bounds.
        for (i, id) in enumerate(var_index_set)
            # If current bounds are within the minimum variable width, skip.
            bounds[i][2] - bounds[i][1] < id[5] && continue

            # Perform optimization-based bound tightening for lower and upper bounds.
            time_elapsed += @elapsed lb = _solve_bound(bt, id, _MOI.MIN_SENSE, bounds[i][1])
            time_elapsed > time_limit && ((terminate = true) && break)
            time_elapsed += @elapsed ub = _solve_bound(bt, id, _MOI.MAX_SENSE, bounds[i][2])
            time_elapsed > time_limit && ((terminate = true) && break)

            # Compute corrected bounds from the optimization-based bounds.
            bounds[i] = _compute_corrected_bounds(lb, ub, bounds[i][1], bounds[i][2], id[5], id[6])
        end

        # Update the iteration counter.
        current_iteration += 1

        # Set the termination variable if max iterations is exceeded.
        current_iteration >= max_iter && (terminate = true)

        # Update the modifications based on the new bounds, then update the original data.
        _update_modifications!(modifications, var_index_set, bounds)
        _IM.update_data!(data, modifications)

        # Get the current average bound widths.
        if use_reduced_network
            avg_widths_final = _get_average_widths(modifications)
        else
            avg_widths_final = _get_average_widths_mn(modifications)
        end

        if maximum(avg_widths_initial - avg_widths_final) < improvement_tol
            terminate = true
        end

        if !terminate
            # Set the initial average widths to the last average widths.
            avg_widths_initial = avg_widths_final

            # Instantiate the new model and set the optimizer.
            bt = instantiate_model(data, model_type, build_type; ext=ext)
            upper_bound_constraint && _constraint_obj_bound(bt, upper_bound)
            relaxed && _relax_model!(bt)
            JuMP.set_optimizer(bt.model, optimizer)
        end
    end

    if use_reduced_network
        _revert_reduced_data!(data)
    end
end


function _check_obbt_options(ub::Float64, rel_gap::Float64, ub_constraint::Bool)
    if ub_constraint && isinf(ub)
        Memento.error(_LOGGER, "[OBBT] The option \"upper_bound_constraint\" cannot be set to true without specifying an upper bound.")
    end

    if !isinf(rel_gap) && isinf(ub)
        Memento.error(_LOGGER, "[OBBT] The option \"rel_gap_tol\" is specified without providing an upper bound.")
    end
end


function _compute_corrected_bounds(lb::Float64, ub::Float64, lb_old::Float64, ub_old::Float64, width::Float64, prec::Int)
    # Convert new bounds to the desired decimal precision.
    lb = max(floor(10.0^prec * lb) * inv(10.0^prec), lb_old)
    ub = min(ceil(10.0^prec * ub) * inv(10.0^prec), ub_old)

    # If the bounds satisfy the minimum width, return them.
    ub - lb >= width && return lb, ub

    # Compute the mean of the bounds.
    mean = 0.5 * (lb + ub)

    # Return bounds satisfying the minimum width.
    if mean - 0.5 * width < lb_old
        return lb_old, ub_old + width
    elseif mean + 0.5 * width > ub_old
        return ub_old - width, ub_old
    else
        return mean - 0.5 * width, mean + 0.5 * width
    end
end


function _solve_bound(wm::AbstractWaterModel, index::Tuple, sense::_MOI.OptimizationSense, start::Float64)
    # Get the variable reference from the index tuple.
    variable = var(wm, index[1], index[3], index[4])

    # Optimize the variable (or affine expression) being tightened.
    JuMP.@objective(wm.model, sense, variable)
    JuMP.optimize!(wm.model)

    # Return an optimized bound or the initial bound that was started with.
    if JuMP.termination_status(wm.model) in [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL]
        # Get the objective value and return the better of the old and new bounds.
        candidate = JuMP.objective_value(wm.model)
        return sense === _MOI.MIN_SENSE ? max(candidate, start) : min(candidate, start)
    else
        message = "[OBBT] Optimization of $(index[1])_$(index[3])[$(index[4])] errored. Adjust tolerances."
        JuMP.termination_status(wm.model) !== _MOI.TIME_LIMIT && Memento.warn(_LOGGER, message)
        return start # Optimization was not successful. Return the starting bound.
    end
end
