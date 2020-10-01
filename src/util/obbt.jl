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


function _get_node_bound_dict(wm::AbstractWaterModel, nw::Int)
    return Dict(string(i) => Dict("h_min"=>JuMP.lower_bound(var(wm, nw, :h, i)),
        "h_max"=>JuMP.upper_bound(var(wm, nw, :h, i))) for i in ids(wm, nw, :node))
end


function _get_edge_bound_dict(wm::AbstractWaterModel, nw::Int, comp::Symbol)
    qp_sym, qn_sym = Symbol("qp_" * string(comp)), Symbol("qn_" * string(comp))
    qp, qn = var(wm, nw, qp_sym), var(wm, nw, qn_sym)

    return Dict(string(a) => Dict("q_min"=>min(0.0, -JuMP.upper_bound(qn[a])),
        "q_max"=>max(0.0, JuMP.upper_bound(qp[a])), "q_min_active"=>0.0)
        for a in ids(wm, nw, comp))
end


function _create_modifications_reduced(wm::AbstractWaterModel)
    data = Dict{String,Any}("per_unit"=>false)
    data["node"] = _get_node_bound_dict(wm, wm.cnw)

    for type in [:pipe, :pump, :regulator, :short_pipe, :valve]
        data[string(type)] = _get_edge_bound_dict(wm, wm.cnw, type)
    end

    return data
end


function _create_modifications_mn(wm::AbstractWaterModel)
    data = Dict{String,Any}("nw"=>Dict{String,Any}(), "per_unit"=>false, "multinetwork"=>true)

    for nw in sort(collect(nw_ids(wm)))
        dnw = data["nw"][string(nw)] = Dict{String,Any}()
        dnw["node"] = _get_node_bound_dict(wm, nw)

        for type in [:pipe, :pump, :regulator, :short_pipe, :valve]
            dnw[string(type)] = _get_edge_bound_dict(wm, nw, type)
        end
    end

    return data
end


function _get_head_index_set(wm::AbstractWaterModel; width::Float64=1.0e-3, prec::Float64=1.0e-4)
    return vcat([[(nw, :node, :h, i, width, prec) for i in ids(wm, nw, :node)]
        for nw in sort(collect(nw_ids(wm)))]...)
end


function _get_flow_index_set(wm::AbstractWaterModel; width::Float64=1.0e-3, prec::Float64=1.0e-4)
    types = [:pipe, :pump, :regulator, :short_pipe, :valve]
    return vcat([vcat([[(nw, type, Symbol("q_" * string(type)), a, width, prec) for a in
        ids(wm, nw, type)] for nw in sort(collect(nw_ids(wm)))]...) for type in types]...)
end


function _get_indicator_index_set(wm::AbstractWaterModel)
    types = [:pump, :regulator, :valve]
    return vcat([vcat([[(nw, type, Symbol("z_" * string(type)), a, 0.0, 0.0) for a in
        ids(wm, nw, type)] for nw in sort(collect(nw_ids(wm)))]...) for type in types]...)
end


function _get_direction_index_set(wm::AbstractWaterModel)
    types = [:pipe, :pump, :regulator, :short_pipe, :valve]
    return vcat([vcat([[(nw, type, Symbol("y_" * string(type)), a, 0.0, 0.0) for a in
        ids(wm, nw, type)] for nw in sort(collect(nw_ids(wm)))]...) for type in types]...)
end


function _get_existing_bounds(data::Dict{String,<:Any}, index::Tuple)
    nw_str, comp_str, index_str = string(index[1]), string(index[2]), string(index[4])
    var_str = occursin("q", string(index[3])) ? "q" : string(index[3])
    min_suffix = comp_str == "pump" ? "_min_active" : "_min"

    if "nw" in keys(data)
        data_index = data["nw"][nw_str][comp_str][index_str]
    else
        data_index = data[comp_str][index_str]
    end

    return data_index[var_str * min_suffix], data_index[var_str * "_max"]
end


function _update_modifications!(data::Dict{String,Any}, original::Dict{String,Any}, var_index_set::Array, bounds::Array)
    # Loop over variables and tighten bounds.
    for (i, id) in enumerate(var_index_set)
        nw_str, comp_str, index_str = string(id[1]), string(id[2]), string(id[4])
        var_str = occursin("q", string(id[3])) ? "q" : string(id[3])
        min_suffix = comp_str == "pump" ? "_min_active" : "_min"

        if _IM.ismultinetwork(data)
            data["nw"][nw_str][comp_str][index_str][var_str * min_suffix] = bounds[i][1]
            data["nw"][nw_str][comp_str][index_str][var_str * "_max"] = bounds[i][2]
        else
            data[comp_str][index_str][var_str * min_suffix] = bounds[i][1]
            data[comp_str][index_str][var_str * "_max"] = bounds[i][2]
        end
    end
end


function _get_average_widths(data::Dict{String,<:Any}, original::Dict{String,<:Any})
    h = sum(x["h_max"] - x["h_min"] for (i, x) in data["node"])
    message = "[OBBT] Average bound widths: h -> $(h * inv(length(data["node"]))), "
    avg_vals = [h * inv(length(data["node"]))]

    for type in ["pipe", "pump", "regulator", "short_pipe", "valve"]
        q_width_sum = 0.0
        min_suffix = type == "pump" ? "_min_active" : "_min"

        for (a, comp) in data[type]
            q_width_sum += comp["q_max"] - comp["q" * min_suffix]
        end

        if length(data[type]) > 0
            q_ave_width = q_width_sum * inv(length(data[type]))
            message *= "q_$(type) -> $(q_ave_width), "
            avg_vals = vcat(avg_vals, q_ave_width)
        end
    end

    Memento.info(_LOGGER, message[1:end-2] * ".")
    return avg_vals
end


function _get_average_widths_mn(data::Dict{String,<:Any}, original::Dict{String,<:Any})
    h = sum(sum(x["h_max"] - x["h_min"] for (i, x) in nw["node"]) for (n, nw) in data["nw"])
    h_length = sum(length(nw["node"]) for (n, nw) in data["nw"])
    message = "[OBBT] Average bound widths: h -> $(h * inv(h_length)), "
    avg_vals = [h * inv(h_length)]

    for type in ["pipe", "pump", "regulator", "short_pipe", "valve"]
        min_suffix = type == "pump" ? "_min_active" : "_min"

        if sum(length(nw[type]) for (n, nw) in data["nw"]) > 0
            q_length = sum(length(nw[type]) for (n, nw) in data["nw"])
            q = sum(sum(x["q_max"] - x["q" * min_suffix] for (i, x) in nw[type]) for (n, nw) in data["nw"])
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


function run_obbt_owf!(data::Dict{String,<:Any}, optimizer; use_reduced_network::Bool = true,
    model_type::Type = CQRDWaterModel, time_limit::Float64 = 3600.0, upper_bound::Float64 =
    Inf, upper_bound_constraint::Bool = false, max_iter::Int = 100, improvement_tol::Float64
    = 1.0e-6, relaxed::Bool = true, precision=1.0e-4, min_width::Float64 = 1.0e-3,
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
    bt = instantiate_model(data, model_type, build_type; ext=ext)
    upper_bound_constraint && _constraint_obj_bound(bt, upper_bound)
    relaxed && _relax_model!(bt) # Relax integrality, if required.

    # Build the dictionary and sets that will store to the network.
    modifications = use_reduced_network ? _create_modifications_reduced(bt) : _create_modifications_mn(bt)
    head_index_set = _get_head_index_set(bt; width=min_width, prec=precision)
    flow_index_set = _get_flow_index_set(bt; width=min_width, prec=precision)
    var_index_set = vcat(head_index_set, flow_index_set)
    bounds = [_get_existing_bounds(modifications, vid) for vid in var_index_set]

    indicator_index_set = _get_indicator_index_set(bt)
    direction_index_set = _get_direction_index_set(bt)

    # Get the initial average bound widths.
    if use_reduced_network
        avg_widths_initial = _get_average_widths(modifications, data)
    else
        avg_widths_initial = _get_average_widths_mn(modifications, data)
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

        # Loop over indicator variables and tighten bounds.
        for (i, id) in enumerate(indicator_index_set)
            # Perform optimization-based bound tightening for lower and upper bounds.
            !JuMP.is_binary(var(bt, id[1], id[3], id[4])) && continue
            time_elapsed += @elapsed lb = _solve_binary_bound(bt, id, _MOI.MIN_SENSE)
            time_elapsed > time_limit && ((terminate = true) && break)
            time_elapsed += @elapsed ub = _solve_binary_bound(bt, id, _MOI.MAX_SENSE)
            time_elapsed > time_limit && ((terminate = true) && break)

            if !_IM.ismultinetwork(data)
                if lb == 1.0 && ub == 1.0
                    modifications[String(id[2])]["$(id[4])"]["force_on"] = true
                elseif lb == 0.0 && ub == 0.0
                    modifications[String(id[2])]["$(id[4])"]["force_off"] = true
                end
            end
        end

        # Loop over direction variables and tighten bounds.
        for (i, id) in enumerate(direction_index_set)
            # Perform optimization-based bound tightening for lower and upper bounds.
            !JuMP.is_binary(var(bt, id[1], id[3], id[4])) && continue
            time_elapsed += @elapsed lb = _solve_binary_bound(bt, id, _MOI.MIN_SENSE)
            time_elapsed > time_limit && ((terminate = true) && break)
            time_elapsed += @elapsed ub = _solve_binary_bound(bt, id, _MOI.MAX_SENSE)
            time_elapsed > time_limit && ((terminate = true) && break)

            if !_IM.ismultinetwork(data)
                if lb == 1.0 && ub == 1.0
                    modifications[String(id[2])]["$(id[4])"]["force_forward"] = true
                elseif lb == 0.0 && ub == 0.0
                    modifications[String(id[2])]["$(id[4])"]["force_reverse"] = true
                end
            end
        end

        # Update the iteration counter.
        current_iteration += 1

        # Set the termination variable if max iterations is exceeded.
        current_iteration >= max_iter && (terminate = true)

        # Update the modifications based on the new bounds, then update the original data.
        _update_modifications!(modifications, data, var_index_set, bounds)
        _IM.update_data!(data, modifications)

        # Get the current average bound widths.
        if use_reduced_network
            avg_widths_final = _get_average_widths(modifications, data)
        else
            avg_widths_final = _get_average_widths_mn(modifications, data)
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


function _check_obbt_options(ub::Float64, ub_constraint::Bool)
    if ub_constraint && isinf(ub)
        Memento.error(_LOGGER, "[OBBT] The option \"upper_bound_constraint\" cannot be set to true without specifying an upper bound.")
    end
end


function _compute_corrected_bounds(lb::Float64, ub::Float64, lb_old::Float64, ub_old::Float64, width::Float64, prec::Float64)
    # Convert new bounds to the desired decimal precision.
    lb = max(floor(inv(prec) * lb) * prec, lb_old)
    ub = min(ceil(inv(prec) * ub) * prec, ub_old)

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


function _solve_binary_bound(wm::AbstractWaterModel, index::Tuple, sense::_MOI.OptimizationSense)
    # Get the variable reference from the index tuple.
    variable = var(wm, index[1], index[3], index[4])

    # Optimize the variable (or affine expression) being tightened.
    JuMP.@objective(wm.model, sense, variable)
    JuMP.optimize!(wm.model)

    # Return an optimized bound or the initial bound that was started with.
    if JuMP.termination_status(wm.model) in [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL]
        # Get the objective value and return the better of the old and new bounds.
        if sense === _MOI.MIN_SENSE
            return JuMP.objective_value(wm.model) >= 0.5 ? 1.0 : 0.0
        elseif sense === _MOI.MAX_SENSE
            return JuMP.objective_value(wm.model) < 0.5 ? 0.0 : 1.0
        end
    else
        message = "[OBBT] Optimization of $(index[1])_$(index[3])[$(index[4])] errored. Adjust tolerances."
        JuMP.termination_status(wm.model) !== _MOI.TIME_LIMIT && Memento.warn(_LOGGER, message)
        return sense == _MOI.MIN_SENSE ? 0.0 : 1.0
    end
end


function _solve_bound(wm::AbstractWaterModel, index::Tuple, sense::_MOI.OptimizationSense, start::Float64)
    # Get the variable reference from the index tuple.
    variable = var(wm, index[1], index[3], index[4])

    # Determine whether or not the associated indicator variable should be fixed.
    fix_status = index[3] == :q_pump ? true : false
  
    if fix_status # Fix the associated indicator variable to one.
        z_var_sym = Symbol(replace(String(index[3]), "q_" => "z_"))
        z_var = var(wm, index[1], z_var_sym, index[4])

        if JuMP.is_binary(z_var)
            JuMP.fix(z_var, 1.0)
        else
            JuMP.fix(z_var, 1.0, force=true)
        end
    end

    # Optimize the variable (or affine expression) being tightened.
    JuMP.@objective(wm.model, sense, variable)
    JuMP.optimize!(wm.model)
    termination_status = JuMP.termination_status(wm.model)

    if termination_status in [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL]
        candidate = JuMP.objective_value(wm.model)
    else
        candidate = nothing
    end

    if fix_status # Reset the indicator variable back to binary.
        z_var_sym = Symbol(replace(String(index[3]), "q_" => "z_"))
        z_var = var(wm, index[1], z_var_sym, index[4])
        JuMP.unfix(z_var)

        if !JuMP.is_binary(z_var)
            JuMP.set_lower_bound(z_var, 0.0)
            JuMP.set_upper_bound(z_var, 1.0)
        end
    end

    # Return an optimized bound or the initial bound that was started with.
    if termination_status in [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL]
        # Get the objective value and return the better of the old and new bounds.
        return sense === _MOI.MIN_SENSE ? max(candidate, start) : min(candidate, start)
    else
        message = "[OBBT] Optimization of $(index[1])_$(index[3])[$(index[4])] errored. Adjust tolerances."
        termination_status !== _MOI.TIME_LIMIT && Memento.warn(_LOGGER, message)
        return start # Optimization was not successful. Return the starting bound.
    end
end
