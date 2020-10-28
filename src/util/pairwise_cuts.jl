mutable struct _PairwiseCutProblem
    sense::_MOI.OptimizationSense
    variable_index_1::_VariableIndex
    variable_index_2::_VariableIndex
    variable_2_fixing_value::Float64
end


function _optimize_bound_problem!(wm::AbstractWaterModel, problem::_PairwiseCutProblem)
    # Get the variables involved in the pairwise problem and fix.
    v_1 = _get_variable_from_index(wm, problem.variable_index_1)
    v_2 = _get_variable_from_index(wm, problem.variable_index_2)
    JuMP.fix(v_2, problem.variable_2_fixing_value) # Fix the second variable.

    # Optimize the first variable (or affine expression).
    JuMP.@objective(wm.model, problem.sense, v_1)
    JuMP.optimize!(wm.model)

    # Return the termination status of the solver.
    return JuMP.termination_status(wm.model)
end


function _get_bound_problem_candidate(wm::AbstractWaterModel, problem::_PairwiseCutProblem)
    if JuMP.termination_status(wm.model) in [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL]
        candidate = JuMP.objective_value(wm.model)

        if problem.sense === _MOI.MIN_SENSE
            return candidate >= 0.5 ? 1.0 : 0.0
        elseif problem.sense === _MOI.MAX_SENSE
            return candidate < 0.5 ? 0.0 : 1.0
        end
    else
        return problem.sense === _MOI.MIN_SENSE ? 0.0 : 1.0
    end
end


function _solve_bound_problem!(wm::AbstractWaterModel, problem::_PairwiseCutProblem)
    termination_status = _optimize_bound_problem!(wm, problem)
    candidate = _get_bound_problem_candidate(wm, problem)
    JuMP.unfix(_get_variable_from_index(wm, problem.variable_index_2))

    if !(termination_status in [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL])
        message = "[OBBT] Optimization of $(problem.variable_index_1) errored."
        termination_status !== _MOI.TIME_LIMIT && Memento.warn(_LOGGER, message)
    end

    return candidate
end


#function _create_pairwise_problems(wm::AbstractWaterModel)
#
#
#    for pipe_id in collect(ids(wm, :pipe))
#        append!(vars, [_VariableIndex(wm.cnw, :pipe, :y_pipe, pipe_id)])
#    end
#end


#function _create_coupled_bound_problems(wm::AbstractWaterModel)
#    vars = _get_indicator_variables(wm)
#    binary_variable_mapping = []
#
#    for id_1 in vars
#        for id_2 in setdiff(vars, [id_1])
#            match_0_min = CoupledBoundProblem(_MOI.MIN_SENSE, id_2, [], [id_1], true)
#            match_0_min_val = _solve_coupled_bound_problem!(wm, match_0_min)
#
#            match_0_max = CoupledBoundProblem(_MOI.MAX_SENSE, id_2, [], [id_1], true)
#            match_0_max_val = _solve_coupled_bound_problem!(wm, match_0_max)
#
#            if match_0_min_val == match_0_max_val
#                append!(binary_variable_mapping, [((0.0, id_1), (match_0_min_val, id_2))])
#            end
#
#            match_1_min = CoupledBoundProblem(_MOI.MIN_SENSE, id_2, [id_1], [], true)
#            match_1_min_val = _solve_coupled_bound_problem!(wm, match_1_min)
#
#            match_1_max = CoupledBoundProblem(_MOI.MAX_SENSE, id_2, [id_1], [], true)
#            match_1_max_val = _solve_coupled_bound_problem!(wm, match_1_max)
#
#            if match_1_min_val == match_1_max_val
#                append!(binary_variable_mapping, [((1.0, id_1), (match_1_min_val, id_2))])
#            end
#        end
#    end
#
#    return binary_variable_mapping
#end

#function _get_binary_variables(wm::AbstractWaterModel)
#    vars = Array{_VariableIndex, 1}()
#
#    for comp_sym in [:pump, :regulator, :valve]
#        comp_str = string(comp_sym)
#        append!(vars, [_VariableIndex(wm.cnw, comp_sym, :z_pump, pump_id)])
#
#    end
#end

#function _get_indicator_variables(wm::AbstractWaterModel)
#    vars = Array{_VariableIndex, 1}()
#
#    for pipe_id in collect(ids(wm, :pipe))
#        append!(vars, [_VariableIndex(wm.cnw, :pipe, :y_pipe, pipe_id)])
#    end
#
#    for pump_id in collect(ids(wm, :pump))
#        append!(vars, [_VariableIndex(wm.cnw, :pump, :y_pump, pump_id)])
#        append!(vars, [_VariableIndex(wm.cnw, :pump, :z_pump, pump_id)])
#    end
#
#    for regulator_id in collect(ids(wm, :regulator))
#        append!(vars, [_VariableIndex(wm.cnw, :regulator, :y_regulator, regulator_id)])
#        append!(vars, [_VariableIndex(wm.cnw, :regulator, :z_regulator, regulator_id)])
#    end
#
#    for short_pipe_id in collect(ids(wm, :short_pipe))
#        append!(vars, [_VariableIndex(wm.cnw, :short_pipe, :y_short_pipe, short_pipe_id)])
#    end
#
#    for valve_id in collect(ids(wm, :valve))
#        append!(vars, [_VariableIndex(wm.cnw, :valve, :y_valve, valve_id)])
#        append!(vars, [_VariableIndex(wm.cnw, :valve, :z_valve, valve_id)])
#    end
#
#    return vars
#end


#mutable struct CoupledBoundProblem
#    sense::_MOI.OptimizationSense
#    variable_to_tighten::_VariableIndex
#    variables_one::Array{_VariableIndex} # Fix to one.
#    variables_zero::Array{_VariableIndex} # Fix to zero.
#    changed::Bool
#end


#function _solve_coupled_bound_problem!(wm::AbstractWaterModel, bound_problem::CoupledBoundProblem)
#    v = _get_variable_from_index(wm, bound_problem.variable_to_tighten)
#    v_one = _get_variable_from_index.(Ref(wm), bound_problem.variables_one)
#    v_zero = _get_variable_from_index.(Ref(wm), bound_problem.variables_zero)
#
#    # Fix binary variables to desired values.
#    JuMP.fix.(v_one, 1.0), JuMP.fix.(v_zero, 0.0)
#
#    # Optimize the variable (or affine expression) being tightened.
#    JuMP.@objective(wm.model, bound_problem.sense, v)
#    JuMP.optimize!(wm.model)
#    termination_status = JuMP.termination_status(wm.model)
#
#    # Store the candidate bound calculated from the optimization.
#    if termination_status in [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL]
#        candidate = JuMP.objective_value(wm.model)
#    end
#
#    # Unfix binary variables that were fixed above.
#    JuMP.unfix.(v_one), JuMP.unfix.(v_zero)
#
#    # Return an optimized bound or the initial bound that was started with.
#    if termination_status in [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL]
#        # Get the objective value and return the better of the old and new bounds.
#        if bound_problem.sense === _MOI.MIN_SENSE
#            return candidate >= 0.5 ? 1.0 : 0.0
#        elseif bound_problem.sense === _MOI.MAX_SENSE
#            return candidate < 0.5 ? 0.0 : 1.0
#        end
#    else
#        message = "[OBBT] Optimization of $(bound_problem.variable_to_tighten) errored. Adjust tolerances."
#        termination_status !== _MOI.TIME_LIMIT && Memento.warn(_LOGGER, message)
#        return bound_problem.sense === _MOI.MIN_SENSE ? 0.0 : 1.0
#    end
#end
#
#
#function _create_coupled_bound_problems(wm::AbstractWaterModel)
#    vars = _get_indicator_variables(wm)
#    binary_variable_mapping = []
#
#    for id_1 in vars
#        for id_2 in setdiff(vars, [id_1])
#            match_0_min = CoupledBoundProblem(_MOI.MIN_SENSE, id_2, [], [id_1], true)
#            match_0_min_val = _solve_coupled_bound_problem!(wm, match_0_min)
#
#            match_0_max = CoupledBoundProblem(_MOI.MAX_SENSE, id_2, [], [id_1], true)
#            match_0_max_val = _solve_coupled_bound_problem!(wm, match_0_max)
#
#            if match_0_min_val == match_0_max_val
#                append!(binary_variable_mapping, [((0.0, id_1), (match_0_min_val, id_2))])
#            end
#
#            match_1_min = CoupledBoundProblem(_MOI.MIN_SENSE, id_2, [id_1], [], true)
#            match_1_min_val = _solve_coupled_bound_problem!(wm, match_1_min)
#
#            match_1_max = CoupledBoundProblem(_MOI.MAX_SENSE, id_2, [id_1], [], true)
#            match_1_max_val = _solve_coupled_bound_problem!(wm, match_1_max)
#
#            if match_1_min_val == match_1_max_val
#                append!(binary_variable_mapping, [((1.0, id_1), (match_1_min_val, id_2))])
#            end
#        end
#    end
#
#    return binary_variable_mapping
#end
#
#
#function compute_binary_cuts(data::Dict{String,<:Any}, optimizer; use_reduced_network::Bool = true,
#    model_type::Type = CQRDWaterModel, time_limit::Float64 = 3600.0, upper_bound::Float64 =
#    Inf, upper_bound_constraint::Bool = false, max_iter::Int = 100, improvement_tol::Float64
#    = 1.0e-6, relaxed::Bool = true, precision = 1.0e-3, min_width::Float64 = 1.0e-2,
#    ext::Dict{Symbol,<:Any} = Dict{Symbol,Any}(:pump_breakpoints=>3), kwargs...)
#    Memento.info(_LOGGER, "[OBBT] Maximum time limit for OBBT set to default value of $(time_limit) seconds.")
#
#    if !_IM.ismultinetwork(data) && use_reduced_network
#        _make_reduced_data!(data)
#        build_type = WaterModels.build_wf
#    elseif _IM.ismultinetwork(data) && !use_reduced_network
#        build_type = WaterModels.build_mn_owf
#    else
#        build_type = WaterModels.build_mn_owf
#        Memento.error(_LOGGER, "[OBBT] Ensure input network is the correct type.")
#    end
#
#    # Check for keyword argument inconsistencies.
#    _check_obbt_options(upper_bound, upper_bound_constraint)
#
#    # Instantiate the bound tightening model and relax integrality, if specified.
#    wm = instantiate_model(data, model_type, build_type; ext=ext)
#    upper_bound_constraint && _constraint_obj_bound(wm, upper_bound)
#    relaxed && _relax_model!(wm) # Relax integrality, if required.
#
#    # Set the optimizer for the bound tightening model.
#    JuMP.set_optimizer(wm.model, optimizer)
#
#    if use_reduced_network
#        _revert_reduced_data!(data)
#    end
#
#    # Get the mapping of binary variables.
#    return _create_coupled_bound_problems(wm)
#end
#
#
#function add_coupled_binary_cuts!(wm::AbstractWaterModel, mappings::Array)
#    num_cuts_added = 0
#
#    for mapping in mappings
#        id_1, id_2 = mapping[1][2], mapping[2][2]
#
#        for nw_1 in sort(collect(nw_ids(wm)))
#            z_1 = var(wm, nw_1, id_1.variable_symbol, id_1.component_index)
#
#            for nw_2 in sort(collect(nw_ids(wm)))
#                z_2 = var(wm, nw_2, id_2.variable_symbol, id_2.component_index)
#
#                if mapping[1][1] == 0.0 && mapping[2][1] == 1.0
#                    JuMP.@constraint(wm.model, z_2 >= 1.0 - z_1)
#                elseif mapping[1][1] == 1.0 && mapping[2][1] == 0.0
#                    JuMP.@constraint(wm.model, z_2 <= 1.0 - z_1)
#                elseif mapping[1][1] == 1.0 && mapping[2][1] == 1.0
#                    JuMP.@constraint(wm.model, z_2 >= z_1)
#                elseif mapping[1][1] == 0.0 && mapping[2][1] == 0.0
#                    JuMP.@constraint(wm.model, z_2 <= z_1)
#                end
#
#                num_cuts_added += 1
#            end
#        end
#    end
#
#    Memento.info(_LOGGER, "[OBBT] Added $(num_cuts_added) coupled binary cuts.")
#end
#
#
#function _sum_node_flows(wm::AbstractWaterModel, nw::Int64, i::Int64)
#    # Collect various indices for node-type components connected to node `i`.
#    reservoirs = ref(wm, nw, :node_reservoir, i) # Reservoirs attached to node `i`.
#    tanks = ref(wm, nw, :node_tank, i) # Tanks attached to node `i`.
#    demands = ref(wm, nw, :node_demand, i) # Demands attached to node `i`.
#
#    # Sum the constant demands required at node `i`.
#    nondispatchable_demands = filter(j -> j in ids(wm, nw, :nondispatchable_demand), demands)
#    fixed_demands = [ref(wm, nw, :nondispatchable_demand, j)["flow_rate"] for j in nondispatchable_demands]
#    net_fixed_demand = length(fixed_demands) > 0 ? sum(fixed_demands) : 0.0
#
#    # Get the indices of dispatchable demands connected to node `i`.
#    dispatchable_demands = filter(j -> j in ids(wm, nw, :dispatchable_demand), demands)
#
#    # Generate the expression for flow.
#    expr = JuMP.AffExpr(net_fixed_demand)
#    expr -= length(tanks) > 0 ? sum(var(wm, nw, :q_tank, i) for i in tanks) : 0.0
#    expr += length(dispatchable_demands) > 0 ? sum(var(wm, nw, :q_demand, i) for i in dispatchable_demands) : 0.0
#    expr -= length(reservoirs) > 0 ? sum(var(wm, nw, :q_reservoir, i) for i in reservoirs) : 0.0
#
#    return expr
#end
#
#
#function _collect_remaining_nws(wm::AbstractWaterModel, n::Int)
#    # Get all network IDs in the multinetwork.
#    network_ids = sort(collect(nw_ids(wm)))
#    n_id = findfirst(x -> x == n, network_ids)
#    return network_ids[n_id:end]
#end
#
#
#function _sum_remaining_pump_flows(wm::AbstractWaterModel, nws::Array{Int64,1})
#    return sum(sum(var(wm, nw, :q_pump)) for nw in nws)
#end
#
#
#function _sum_remaining_demands(wm::AbstractWaterModel, nws::Array{Int64,1})
#    expr = JuMP.AffExpr(0.0)
#
#    for nw in nws
#        for (i, node) in ref(wm, nw, :node)
#            # Sum the constant demands required at node `i`.
#            demands = ref(wm, nw, :node_demand, i) # Demands attached to node `i`.
#            nondispatchable_demands = filter(j -> j in ids(wm, nw, :nondispatchable_demand), demands)
#            fixed_demands = [ref(wm, nw, :nondispatchable_demand, j)["flow_rate"] for j in nondispatchable_demands]
#            net_fixed_demand = length(fixed_demands) > 0 ? sum(fixed_demands) : 0.0
#
#            # Get the indices of dispatchable demands connected to node `i`.
#            dispatchable_demands = filter(j -> j in ids(wm, nw, :dispatchable_demand), demands)
#
#            # Generate the expression for flow.
#            expr += length(dispatchable_demands) > 0 ? sum(var(wm, nw, :q_demand, i) for i in dispatchable_demands) : 0.0
#            expr += length(fixed_demands) > 0 ? sum(fixed_demands) : 0.0
#        end
#    end
#
#    return expr
#end
#
#
#function _add_capacity_cuts!(wm::AbstractWaterModel)
#    # Get all network IDs in the multinetwork.
#    network_ids = sort(collect(nw_ids(wm)))
#
#    # Start with the first network, representing the initial time step.
#    n_1 = network_ids[1]
#
#    for (n, network) in nws(wm)
#        nws_remaining = _collect_remaining_nws(wm, n)
#        time_step = ref(wm, n, :time_step)
#        demand_volume_remaining = _sum_remaining_demands(wm, nws_remaining) * time_step
#        pump_volume_remaining = _sum_remaining_pump_flows(wm, nws_remaining) * time_step
#        tank_volume_sum = sum(var(wm, n_1, :V)) - sum(var(wm, n, :V))
#        c = JuMP.@constraint(wm.model, tank_volume_sum + demand_volume_remaining <= pump_volume_remaining)
#    end
#end
