mutable struct _PairwiseProblem
    sense::JuMP.OptimizationSense
    variable_index_1::_VariableIndex
    variable_index_2::_VariableIndex
    variable_2_fixing_value::Float64
    precision::Float64
end


mutable struct _PairwiseProblemSet
    problems::Array{_PairwiseProblem, 1}
end


mutable struct _PairwiseCut
    coefficient_1::Float64
    variable_index_1::_VariableIndex
    coefficient_2::Float64
    variable_index_2::_VariableIndex
    constant::Float64
end


function _optimize_bound_problem!(wm::AbstractWaterModel, problem::_PairwiseProblem)
    # Get the variables involved in the pairwise problem and fix.
    v_1 = _get_variable_from_index(wm, problem.variable_index_1)
    v_2 = _get_variable_from_index(wm, problem.variable_index_2)

    if !JuMP.is_fixed(v_2)
        # Fix the second variable in the pairwise problem (always binary).
        JuMP.fix(v_2, problem.variable_2_fixing_value)
    end

    # Optimize the first variable (or affine expression).
    JuMP.@objective(wm.model, problem.sense, v_1)
    JuMP.optimize!(wm.model)

    # Return the termination status of the solver.
    return JuMP.termination_status(wm.model)
end


function _get_bound_problem_candidate(wm::AbstractWaterModel, problem::_PairwiseProblem)
    v_1 = _get_variable_from_index(wm, problem.variable_index_1)
    v_1_lb = _get_lower_bound_from_index(wm, problem.variable_index_1)
    v_1_ub = _get_upper_bound_from_index(wm, problem.variable_index_1)
    termination_status = JuMP.termination_status(wm.model)

    # Otherwise, store the original bound.
    if problem.sense === JuMP.MIN_SENSE
        candidate = v_1_lb
    elseif problem.sense === JuMP.MAX_SENSE
        candidate = v_1_ub
    end

    if !(termination_status in [JuMP.INFEASIBLE, JuMP.INFEASIBLE_OR_UNBOUNDED])
        try
            # Store the candidate value calculated from the optimization.
            candidate = JuMP.objective_bound(wm.model)
        catch
            # Otherwise, store the original bound.
            if problem.sense === JuMP.MIN_SENSE
                candidate = v_1_lb
            elseif problem.sense === JuMP.MAX_SENSE
                candidate = v_1_ub
            end
        end

        if isa(v_1, JuMP.VariableRef) && JuMP.is_binary(v_1)
            if problem.sense === JuMP.MIN_SENSE
                return candidate >= 0.5 ? 1.0 : 0.0
            elseif problem.sense === JuMP.MAX_SENSE
                return candidate < 0.5 ? 0.0 : 1.0
            end
        else
            prec = problem.precision

            if problem.sense === JuMP.MIN_SENSE
                scaled_candidate = floor(inv(prec) * candidate) * prec
                return max(scaled_candidate, v_1_lb)
            elseif problem.sense === JuMP.MAX_SENSE
                scaled_candidate = ceil(inv(prec) * candidate) * prec
                return min(scaled_candidate, v_1_ub)
            end
        end
    else
        if isa(v_1, JuMP.VariableRef) && JuMP.is_binary(v_1)
            return problem.sense === JuMP.MIN_SENSE ? 0.0 : 1.0
        else
            return problem.sense === JuMP.MIN_SENSE ? v_1_lb : v_1_ub
        end
    end
end


function _unfix_bound_problem_variable!(wm::AbstractWaterModel, problem::_PairwiseProblem)
    v = _get_variable_from_index(wm, problem.variable_index_2)

    if JuMP.is_fixed(v)
        JuMP.unfix(v)
    end
end


function _solve_bound_problem!(wm::AbstractWaterModel, problem::_PairwiseProblem)
    termination_status = _optimize_bound_problem!(wm, problem)
    candidate = _get_bound_problem_candidate(wm, problem)
    _unfix_bound_problem_variable!(wm, problem)
    return candidate
end


function _get_pairwise_problem_set(variable_index_1::_VariableIndex, variable_index_2::_VariableIndex, fixing_value::Float64, precision::Float64)
    problem_1 = _PairwiseProblem(JuMP.MIN_SENSE, variable_index_1, variable_index_2, fixing_value, precision)
    problem_2 = _PairwiseProblem(JuMP.MAX_SENSE, variable_index_1, variable_index_2, fixing_value, precision)
    return _PairwiseProblemSet([problem_1, problem_2])
end


function _get_pairwise_problem_sets(wm::AbstractWaterModel; nw::Int = nw_id_default)
    problem_sets = Array{_PairwiseProblemSet, 1}()
    binary_variable_indices = _get_binary_variable_indices(wm; nw = nw)

    flow_variable_indices = _get_flow_variable_indices(wm; nw = nw)
    head_variable_indices = _get_head_variable_indices(wm; nw = nw)
    continuous_variable_indices = vcat(flow_variable_indices, head_variable_indices)

    for variable_index_1 in vcat(binary_variable_indices, continuous_variable_indices)
        for variable_index_2 in setdiff(binary_variable_indices, [variable_index_1])
            problem_set_0 = _get_pairwise_problem_set(variable_index_1, variable_index_2, 0.0, 1.0e-5)
            problem_set_1 = _get_pairwise_problem_set(variable_index_1, variable_index_2, 1.0, 1.0e-5)
            append!(problem_sets, [problem_set_0, problem_set_1])
        end
    end

    return problem_sets
end


function _get_binary_continuous_cut(problem::_PairwiseProblem, value::Float64, v_1_lb::Float64, v_1_ub::Float64)
    vid_1, vid_2 = problem.variable_index_1, problem.variable_index_2

    if problem.sense === JuMP.MIN_SENSE && problem.variable_2_fixing_value == 1.0
        # value * v_2 + v_1_lb * (1.0 - v_2) <= v_1  -->  -v_1 + (value - v_1_lb) * v_2 + v_1_lb <= 0
        return _PairwiseCut(-1.0, vid_1, value - v_1_lb, vid_2, v_1_lb)
    elseif problem.sense === JuMP.MIN_SENSE && problem.variable_2_fixing_value == 0.0
        # value * (1.0 - v_2) + v_1_lb * v_2 <= v_1  -->  -v_1 + (v_1_lb - value) * v_2 + value <= 0
        return _PairwiseCut(-1.0, vid_1, v_1_lb - value, vid_2, value)
    elseif problem.sense === JuMP.MAX_SENSE && problem.variable_2_fixing_value == 1.0
        # v_1 <= value * v_2 + v_1_ub * (1.0 - v_2)  -->  v_1 - (value - v_1_ub) * v_2 - v_1_ub <= 0
        return _PairwiseCut(1.0, vid_1, v_1_ub - value, vid_2, -v_1_ub)
    elseif problem.sense === JuMP.MAX_SENSE && problem.variable_2_fixing_value == 0.0
        # v_1 <= value * (1.0 - v_2) + v_1_ub * v_2  -->  v_1 + (value - v_1_ub) * v_2 - value <= 0
        return _PairwiseCut(1.0, vid_1, value - v_1_ub, vid_2, -value)
    end
end


function _get_binary_pairwise_cut(problem_set::_PairwiseProblemSet, variable_1_value::Float64)
    vid_1 = problem_set.problems[1].variable_index_1
    vid_2 = problem_set.problems[1].variable_index_2

    if variable_1_value == 1.0 && problem_set.problems[1].variable_2_fixing_value == 1.0
        return _PairwiseCut(-1.0, vid_1, 1.0, vid_2, 0.0) # -v_1 + v_2 + 0.0 <= 0.0
    elseif variable_1_value == 1.0 && problem_set.problems[1].variable_2_fixing_value == 0.0
        return _PairwiseCut(-1.0, vid_1, -1.0, vid_2, 1.0) # -v_1 - v_2 + 1.0 <= 0.0
    elseif variable_1_value == 0.0 && problem_set.problems[1].variable_2_fixing_value == 1.0
        return _PairwiseCut(1.0, vid_1, 1.0, vid_2, -1.0) # v_1 + v_2 - 1.0 <= 0.0
    elseif variable_1_value == 0.0 && problem_set.problems[1].variable_2_fixing_value == 0.0
        return _PairwiseCut(1.0, vid_1, -1.0, vid_2, 0.0) # v_1 - v_2 + 0.0 <= 0.0
    end
end


function _compute_pairwise_cuts!(wm::AbstractWaterModel, problem_sets::Array{_PairwiseProblemSet, 1})
    cuts = Vector{_PairwiseCut}([])

    for problem_set in problem_sets
        # Solve all the problems in the _PairwiseProblemSet.
        candidates = _solve_bound_problem!.(Ref(wm), problem_set.problems)
        v_1_vars = [_get_variable_from_index(wm, x.variable_index_1) for x in problem_set.problems]
        v_1_lb = _get_lower_bound_from_index(wm, problem_set.problems[1].variable_index_1)
        v_1_ub = _get_upper_bound_from_index(wm, problem_set.problems[1].variable_index_1)

        if all(x -> isa(x, JuMP.VariableRef) && JuMP.is_binary(x), v_1_vars) && all(x -> x == candidates[1], candidates)
            # If all the candidate solutions are equal, a cut can be added.
            append!(cuts, [_get_binary_pairwise_cut(problem_set, candidates[1])])
        elseif !(all(x -> isa(x, JuMP.VariableRef) && JuMP.is_binary(x), v_1_vars))
            # Add a cut relating a binary variable to a continuous variable.
            for (i, problem) in enumerate(problem_set.problems)
                if (problem.sense === JuMP.MIN_SENSE && candidates[i] > v_1_lb)
                    push!(cuts, _get_binary_continuous_cut(problem, candidates[i], v_1_lb, v_1_ub))
                elseif (problem.sense === JuMP.MAX_SENSE && candidates[i] < v_1_ub)
                    push!(cuts, _get_binary_continuous_cut(problem, candidates[i], v_1_lb, v_1_ub))
                end
            end
        end
    end

    return cuts
end


function _add_pairwise_cuts!(wm::AbstractWaterModel, cuts::Vector{_PairwiseCut})
    for cut in cuts
        has_v_1 = _model_has_variable_index(wm, cut.variable_index_1)
        has_v_2 = _model_has_variable_index(wm, cut.variable_index_2)

        if has_v_1 && has_v_2
            v_1 = _get_variable_from_index(wm, cut.variable_index_1)
            v_2 = _get_variable_from_index(wm, cut.variable_index_2)
            expr = cut.coefficient_1 * v_1 + cut.coefficient_2 * v_2 + cut.constant
            JuMP.@constraint(wm.model, expr <= 0.0)
        end
    end
end


function _add_pairwise_cuts!(wm::AbstractWaterModel; nw::Int = nw_id_default)
    problem_sets = _get_pairwise_problem_sets(wm; nw=nw)
    cuts = _compute_pairwise_cuts!(wm, problem_sets)
    _add_pairwise_cuts!(wm, cuts)
end


function _read_pairwise_cuts(path::String)
    cuts_array = Vector{_PairwiseCut}([])
    
    for entry in JSON.parsefile(path)
        vid_1_network_index = Int(entry["variable_index_1"]["network_index"])
        vid_1_component_type = Symbol(entry["variable_index_1"]["component_type"])
        vid_1_variable_symbol = Symbol(entry["variable_index_1"]["variable_symbol"])
        vid_1_component_index = Int(entry["variable_index_1"]["component_index"])
        vid_1 = _VariableIndex(vid_1_network_index, vid_1_component_type,
            vid_1_variable_symbol, vid_1_component_index)

        vid_2_network_index = Int(entry["variable_index_2"]["network_index"])
        vid_2_component_type = Symbol(entry["variable_index_2"]["component_type"])
        vid_2_variable_symbol = Symbol(entry["variable_index_2"]["variable_symbol"])
        vid_2_component_index = Int(entry["variable_index_2"]["component_index"])
        vid_2 = _VariableIndex(vid_2_network_index, vid_2_component_type,
            vid_2_variable_symbol, vid_2_component_index)

        coefficient_1 = entry["coefficient_1"]
        coefficient_2 = entry["coefficient_2"]
        constant = entry["constant"]
        
        push!(cuts_array, _PairwiseCut(coefficient_1, vid_1,
            coefficient_2, vid_2, constant))
    end

    return cuts_array
end