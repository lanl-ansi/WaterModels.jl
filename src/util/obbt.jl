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
* `min_bound_width`: domain beyond which bound tightening is not performed.
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

function _get_obbt_node_var_lb(wm::AbstractWaterModel, v::Symbol)
    network_ids = sort(collect(nw_ids(wm)))
    return Dict(nw => Dict(i => JuMP.lower_bound(var(wm, nw, v, i))
        for i in ids(wm, nw, :node)) for nw in network_ids)
end

function _get_obbt_node_var_ub(wm::AbstractWaterModel, v::Symbol)
    network_ids = sort(collect(nw_ids(wm)))
    return Dict(nw => Dict(i => JuMP.upper_bound(var(wm, nw, v, i))
        for i in ids(wm, nw, :node)) for nw in network_ids)
end

function run_obbt_owf!(data::Dict{String,<:Any}, optimizer;
    model_type::Type = MICPRWaterModel,
    max_iter::Int = 100,
    time_limit::Float64 = 3600.0,
    upper_bound::Float64 = Inf,
    upper_bound_constraint::Bool = false,
    rel_gap_tol::Float64 = Inf,
    min_bound_width::Float64 = 1.0e-2,
    improvement_tol::Float64 = 1.0e-3,
    precision::Int = 4,
    termination::Symbol = :avg,
    kwargs...)

    Memento.info(_LOGGER, "Maximum number of OBBT iterations set to default value of $(max_iter).")
    Memento.info(_LOGGER, "Maximum time limit for OBBT set to default value of $(time_limit) seconds.")
    model_relaxation = instantiate_model(data, model_type, WaterModels.build_mn_owf)

    # Check for model_type compatability with OBBT.
    _check_variables(model_relaxation)

    # Check for other keyword argument consistencies.
    _check_obbt_options(upper_bound, rel_gap_tol, upper_bound_constraint)

    # Check termination norm criteria for OBBT.
    if termination != :avg && termination != :max
        Memento.error(_LOGGER, "OBBT termination criteria can only be :max or :avg.")
    end

    # Further relax all binary variables in the model.
    bin_vars = filter(v -> JuMP.is_binary(v), JuMP.all_variables(model_relaxation.model))
    JuMP.unset_binary.(bin_vars) # Make all binary variables free and continuous.
    JuMP.set_lower_bound.(bin_vars, 0.0) # Lower-bound the relaxed binary variables.
    JuMP.set_upper_bound.(bin_vars, 1.0) # Upper-bound the relaxed binary variables.

    # Declare possible pass statuses.
    pass_statuses = [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL]

    # Compute initial relative gap between relaxation objective and upper_bound.
    result_relaxation = optimize_model!(model_relaxation, optimizer=optimizer)
    current_relaxation_objective = result_relaxation["objective"]

    if upper_bound < current_relaxation_objective
        Memento.error(_LOGGER, "The upper bound provided to OBBT is not a valid OWF upper bound.")
    end

    if !(result_relaxation["termination_status"] in pass_statuses)
        termination_status = result_relaxation["termination_status"]
        Memento.warn(_LOGGER, "Initial relaxation solve status is $(termination_status).")

        if termination_status == :SubOptimal
            Memento.warn(_LOGGER, "Continuing with the bound tightening algorithm.")
        end
    end

    current_rel_gap = Inf

    if !isinf(upper_bound)
        current_rel_gap = (upper_bound - current_relaxation_objective) * inv(upper_bound)
        Memento.info(_LOGGER, "Initial relaxation gap is $(current_rel_gap).")
    end

    model_bt = instantiate_model(data, model_type, WaterModels.build_mn_owf)
    upper_bound_constraint && _constraint_obj_bound(model_bt, upper_bound)

    # Further relax all binary variables in the bound tightening model.
    bin_vars = filter(v -> JuMP.is_binary(v), JuMP.all_variables(model_bt.model))
    JuMP.unset_binary.(bin_vars) # Make all binary variables free and continuous.
    JuMP.set_lower_bound.(bin_vars, 0.0) # Lower-bound the relaxed binary variables.
    JuMP.set_upper_bound.(bin_vars, 1.0) # Upper-bound the relaxed binary variables.

    # Initialize the statistics dictionary.
    stats = Dict{String,Any}()
    stats["model_type"] = model_type
    stats["initial_relaxation_objective"] = current_relaxation_objective
    stats["initial_rel_gap_from_ub"] = current_rel_gap
    stats["upper_bound"] = upper_bound

    h_lb = _get_obbt_node_var_lb(model_bt, :h)
    h_ub = _get_obbt_node_var_ub(model_bt, :h)
    stats["h_range_init"], stats["avg_h_range_init"] = 0.0, 0.0

    for nw in sort(collect(nw_ids(model_bt)))
        node_ids = ids(model_bt, nw, :node)
        h_range_nw = sum(h_ub[nw][i] - h_lb[nw][i] for i in node_ids)
        stats["h_range_init"] += h_range_nw
        stats["avg_h_range_init"] += h_range_nw * inv(length(node_ids))
    end

    # Initialize algorithm termination metadata.
    h_range_final, total_h_reduction = 0.0, Inf
    max_h_reduction, avg_h_reduction = Inf, Inf
    final_relaxation_objective = NaN
    current_iteration, time_elapsed = 0, 0.0
    parallel_time_elapsed, terminate = 0.0, false

    # Set whether or not the algorithm should immediately terminate.
    if termination == :avg
        terminate = avg_h_reduction <= improvement_tol && avg_q_reduction <= improvement_tol
    elseif termination == :max
        terminate = max_h_reduction <= improvement_tol && max_q_reduction <= improvement_tol
    end

    while !terminate # Algorithmic loop.
        # Set important metadata describing the current round.
        iter_start_time, max_h_iteration_time = time(), 0.0
        total_h_reduction, avg_h_reduction, max_h_reduction = 0.0, 0.0, 0.0

        # Loop over all subnetworks in the multinetwork.
        for nw in sort(collect(nw_ids(model_bt)))

            # Loop over all nodes in the network.
            for i in ids(model_bt, nw, :node)
                println(nw, " ", i)

                # If the lower and upper bounds are already close, skip.
                if h_ub[nw][i] - h_lb[nw][i] < min_bound_width
                    continue
                end

                # Store metadata for variable tightening.
                lb, ub, start_time = NaN, NaN, time()

                # Minimize the variable whose bounds are being tightened.
                JuMP.@objective(model_bt.model, _MOI.MIN_SENSE, var(model_bt, nw, :h, i))
                result_bt = optimize_model!(model_bt, optimizer=optimizer)

                # Store the new lower bound or print a warning and continue.
                if result_bt["termination_status"] in pass_statuses
                    # Compute the candidate lower bound at the specified precision.
                    val = JuMP.objective_value(model_bt.model)
                    lb_new = floor(10.0^precision * val) * inv(10.0^precision)
                    lb_new > h_lb[nw][i] && (lb = lb_new) # Store the new lower bound.
                else
                    Memento.warn(_LOGGER, "BT minimization for h[$(i)] errored. Adjust tolerances.")
                    continue
                end

                # Maximize the variable whose bounds are being tightened.
                JuMP.@objective(model_bt.model, _MOI.MAX_SENSE, var(model_bt, nw, :h, i))
                result_bt = optimize_model!(model_bt, optimizer=optimizer)

                # Store the new lower bound or print a warning and continue.
                if result_bt["termination_status"] in pass_statuses
                    # Compute the candidate lower bound at the specified precision.
                    val = JuMP.objective_value(model_bt.model)
                    ub_new = ceil(10.0^precision * val) * inv(10.0^precision)
                    ub_new < h_ub[nw][i] && (ub = ub_new) # Store the new lower bound.
                else
                    Memento.warn(_LOGGER, "BT maximization for h[$(i)] errored. Adjust tolerances.")
                    continue
                end

                # Compute the time required for the variable's bound tightening.
                end_time = time() - start_time
                max_h_iteration_time = max(end_time, max_h_iteration_time)

                # Perform sanity checks on the new bounds.
                if lb > ub # Check if the lower bound exceeds the upper bound.
                    Memento.warn(_LOGGER, "BT lb > ub for h[$(i)]. Adjust tolerances.")
                    continue
                end

                # Revert to old bounds if new bounds are nonsensical or NaN.
                (!isnan(lb) && lb > h_ub[i]) && (lb = h_lb[i])
                (!isnan(ub) && ub < h_lb[i]) && (ub = h_ub[i])
                isnan(lb) && (lb = h_lb[i])
                isnan(ub) && (ub = h_ub[i])
            end
        end
    end

    #        # vm bound-reduction computation
    #        h_reduction = 0.0
    #        if (ub - lb >= min_bound_width)
    #            h_reduction = (h_ub[node] - h_lb[node]) - (ub - lb)
    #            h_lb[node] = lb
    #            h_ub[node] = ub
    #        else
    #            mean = 0.5 * (ub + lb)
    #            if (mean - 0.5 * min_bound_width < h_lb[node])
    #                lb = h_lb[node]
    #                ub = h_lb[node] + min_bound_width
    #            elseif (mean + 0.5 * min_bound_width > h_ub[node])
    #                ub = h_ub[node]
    #                lb = h_ub[node] - min_bound_width
    #            else
    #                lb = mean - 0.5 * min_bound_width
    #                ub = mean + 0.5 * min_bound_width
    #            end
    #            h_reduction = (h_ub[node] - h_lb[node]) - (ub - lb)
    #            h_lb[node] = lb
    #            h_ub[node] = ub
    #        end

    #        total_h_reduction += (h_reduction)
    #        max_h_reduction = max(h_reduction, max_h_reduction)
    #    end
    #    avg_h_reduction = total_h_reduction/length(nodes)

    #    h_range_final = sum([h_ub[node] - h_lb[node] for node in nodes])

    #    q_range_final = sum([q_ub[bp] - q_lb[bp] for bp in pipes])

    #    parallel_time_elapsed += max(max_h_iteration_time, max_q_iteration_time)

    #    time_elapsed += (time() - iter_start_time)

    #    # populate the modifications, update the data, and rebuild the bound tightening model
    #    modifications = _create_modifications(model_bt, h_lb, h_ub, q_lb, q_ub)
    #    WaterModels.update_data!(data, modifications)
    #    model_bt = instantiate_model(data, model_type, WaterModels.build_mn_owf)
    #    (upper_bound_constraint) && (_constraint_obj_bound(model_bt, upper_bound))
    #    vm = var(model_bt, :vm)
    #    q = var(model_bt, :q)

    #    # run the qc relaxation for the updated bounds
    #    result_relaxation = run_owf(data, model_type::Type, optimizer)

    #    if result_relaxation["termination_status"] in pass_statuses
    #        current_rel_gap = (upper_bound - result_relaxation["objective"])/upper_bound
    #        final_relaxation_objective = result_relaxation["objective"]
    #    else
    #        Memento.warn(_LOGGER, "relaxation solve failed in iteration $(current_iteration+1)")
    #        Memento.warn(_LOGGER, "using the previous iteration's gap to check relative gap stopping criteria")
    #    end

    #    Memento.info(_LOGGER, "iteration $(current_iteration+1), vm range: $h_range_final, q range: $q_range_final, relaxation obj: $final_relaxation_objective")

    #    # Update the termination flag.
    #    if termination == :avg
    #        terminate = avg_h_reduction <= improvement_tol && avg_q_reduction <= improvement_tol
    #    else
    #        terminate = max_h_reduction <= improvement_tol && max_q_reduction <= improvement_tol
    #    end

    #    # iteration counter update
    #    current_iteration += 1
    #    # check all the stopping criteria
    #    (current_iteration >= max_iter) && (Memento.info(_LOGGER, "maximum iteration limit reached"); break)
    #    (time_elapsed > time_limit) && (Memento.info(_LOGGER, "maximum time limit reached"); break)
    #    if (!isinf(rel_gap_tol)) && (current_rel_gap < rel_gap_tol)
    #        Memento.info(_LOGGER, "relative optimality gap < $rel_gap_tol")
    #        break
    #    end

    #end

    #pipes_vad_same_sign_count = 0
    #for (key, pipe) in data["pipe"]
    #    is_same_sign = (pipe["q_max"] >=0 && pipe["q_min"] >= 0) || (pipe["q_max"] <=0 && pipe["q_min"] <= 0)
    #    (is_same_sign) && (pipes_vad_same_sign_count += 1)
    #end

    #stats["final_relaxation_objective"] = final_relaxation_objective
    #stats["final_rel_gap_from_ub"] = isnan(upper_bound) ? Inf : current_rel_gap
    #stats["h_range_final"] = h_range_final
    #stats["avg_h_range_final"] = h_range_final * inv(length(nodes))

    #stats["q_range_final"] = q_range_final
    #stats["avg_q_range_final"] = q_range_final * inv(length(pipes))

    #stats["run_time"] = time_elapsed
    #stats["iteration_count"] = current_iteration
    #stats["sim_parallel_run_time"] = parallel_time_elapsed

    #stats["vad_sign_determined"] = pipes_vad_same_sign_count

    return data, stats
end


function _check_variables(wm::AbstractWaterModel)
    try
        h = var(wm, :h)
    catch err
        (isa(error, KeyError)) && (Memento.error(_LOGGER, "OBBT is not supported for models without head variables."))
    end
end


function _check_obbt_options(ub::Float64, rel_gap::Float64, ub_constraint::Bool)
    if ub_constraint && isinf(ub)
        Memento.error(_LOGGER, "The option \"upper_bound_constraint\" cannot be set to true without specifying an upper bound.")
    end

    if !isinf(rel_gap) && isinf(ub)
        Memento.error(_LOGGER, "The option \"rel_gap_tol\" is specified without providing an upper bound.")
    end
end


function _constraint_obj_bound(wm::AbstractWaterModel, bound)
end


function _create_modifications(wm::AbstractWaterModel,
    h_lb::Dict{Any,Float64}, h_ub::Dict{Any,Float64},
    q_lb::Dict{Any,Float64}, q_ub::Dict{Any,Float64})

    modifications = Dict{String, Any}()
    modifications["per_unit"] = false
    modifications["node"] = Dict{String, Any}()
    modifications["pipe"] = Dict{String, Any}()

    for node in ids(wm, :node)
        index = string(ref(wm, :node, node, "index"))
        modifications["node"][index] = Dict("h_min" => h_lb[node], "h_max" => h_ub[node])
    end

    return modifications
end
