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


function _create_modifications(wm::AbstractWaterModel)
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


function _get_head_index_set(wm::AbstractWaterModel; width::Float64=0.1, prec::Int=4)
    return vcat([[(nw, :node, :h, i, width, prec) for i in ids(wm, nw, :node)]
        for nw in sort(collect(nw_ids(wm)))]...)
end


function _get_flow_index_set(wm::AbstractWaterModel; width::Float64=1.0e-5, prec::Int=7)
    types = [:pipe, :shutoff_valve, :check_valve, :pressure_reducing_valve, :pump]
    return vcat([vcat([[(nw, type, Symbol("q_" * string(type)), a, width, prec) for a in
        ids(wm, nw, type)] for nw in sort(collect(nw_ids(wm)))]...) for type in types]...)
end


function _get_existing_bounds(data::Dict{String,<:Any}, index::Tuple)
    nw_str, comp_str, index_str = string(index[1]), string(index[2]), string(index[4])
    comp_str = comp_str in ["check_valve", "shutoff_valve"] ? "pipe" : comp_str
    var_str = occursin("q", string(index[3])) ? "q" : string(index[3])
    data_index = data["nw"][nw_str][comp_str][index_str]
    return data_index[var_str * "_min"], data_index[var_str * "_max"]
end


function _update_modifications!(data::Dict{String,Any}, var_index_set::Array, bounds::Array)
    # Loop over variables and tighten bounds.
    for (i, id) in enumerate(var_index_set)
        nw_str, comp_str, index_str = string(id[1]), string(id[2]), string(id[4])
        var_str = occursin("q", string(id[3])) ? "q" : string(id[3])
        comp_c = comp_str in ["shutoff_valve", "check_valve"] ? "pipe" : comp_str
        data["nw"][nw_str][comp_c][index_str][var_str * "_min"] = bounds[i][1]
        data["nw"][nw_str][comp_c][index_str][var_str * "_max"] = bounds[i][2]
    end
end


function _print_average_widths(data::Dict{String,<:Any})
    h = sum(sum(x["h_max"] - x["h_min"] for (i, x) in nw["node"]) for (n, nw) in data["nw"])
    h_length = sum(length(nw["node"]) for (n, nw) in data["nw"])
    message = "[OBBT] Average bound widths: h -> $(h * inv(h_length)), "

    for type in ["pipe", "pressure_reducing_valve", "pump"]
        if sum(length(nw[type]) for (n, nw) in data["nw"]) > 0
            q_length = sum(length(nw[type]) for (n, nw) in data["nw"])
            q = sum(sum(x["q_max"] - x["q_min"] for (i, x) in nw[type]) for (n, nw) in data["nw"])
            message *= "q_$(type) -> $(q * inv(q_length)), "
        end
    end

    Memento.info(_LOGGER, message[1:end-2] * ".")
end


function run_obbt_owf!(data::Dict{String,<:Any}, optimizer;
    model_type::Type = MICPRWaterModel, time_limit::Float64 = 3600.0, upper_bound::Float64 =
    Inf, upper_bound_constraint::Bool = false, rel_gap_tol::Float64 = Inf,
    improvement_tol::Float64 = 1.0e-3, termination::Symbol = :avg, relaxed::Bool = true,
    ext::Dict{Symbol,<:Any} = Dict{Symbol,Any}(:pump_breakpoints=>5), kwargs...)
    # Print a message with relevant algorithm limit information.
    Memento.info(_LOGGER, "[OBBT] Maximum time limit for OBBT set to default value of $(time_limit) seconds.")

    # Check for keyword argument inconsistencies.
    _check_obbt_options(upper_bound, rel_gap_tol, upper_bound_constraint)

    # Instantiate the bound tightening model and relax integrality, if specified.
    bt = instantiate_model(data, model_type, WaterModels.build_mn_owf; ext=ext)
    upper_bound_constraint && _constraint_obj_bound(bt, upper_bound)
    relaxed && _relax_model!(bt)

    # Build the dictionary and sets that will store to the network.
    modifications = _create_modifications(bt)
    var_index_set = vcat(_get_head_index_set(bt), _get_flow_index_set(bt))
    bounds = [_get_existing_bounds(modifications, vid) for vid in var_index_set]

    # Print the initial average bound widths.
    _print_average_widths(modifications)

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

        # Update the modifications based on the new bounds, then update the original data.
        _update_modifications!(modifications, var_index_set, bounds)
        _IM.update_data!(data, modifications)

        # Print the current average bound widths.
        _print_average_widths(data)

        # Instantiate the new model and set the optimizer.
        bt = instantiate_model(data, model_type, build_mn_owf; ext=ext)
        upper_bound_constraint && _constraint_obj_bound(bt, upper_bound)
        relaxed && _relax_model!(bt)
        JuMP.set_optimizer(bt.model, optimizer)
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
        message = "[OBBT] Optimization of $(index_1)_$(index[3])[$(index[4])] errored. Adjust tolerances."
        JuMP.termination_status(wm.model) !== _MOI.TIME_LIMIT && Memento.warn(_LOGGER, message)
        return start # Optimization was not successful. Return the starting bound.
    end
end
