function _get_binary_cut_mapping(data::Dict{String,<:Any}, optimizer; model_type::Type = CQRDWaterModel, ext::Dict{Symbol,<:Any} = Dict{Symbol,Any}())
    _make_reduced_data!(data)
    wm = instantiate_model(data, model_type, build_wf; ext=ext)
    indicator_index_set = _get_indicator_index_set(wm)
    JuMP.set_optimizer(wm.model, optimizer)
    binary_index_set = vcat(_get_indicator_index_set(wm), _get_direction_index_set(wm))
    binary_variable_mapping = []

    for id_1 in binary_index_set
        for id_2 in setdiff(binary_index_set, [id_1])
            match_0_min = _solve_coupled_binary(wm, id_1, id_2, 0.0, _MOI.MIN_SENSE)
            match_0_max = _solve_coupled_binary(wm, id_1, id_2, 0.0, _MOI.MAX_SENSE)
            match_1_min = _solve_coupled_binary(wm, id_1, id_2, 1.0, _MOI.MIN_SENSE)
            match_1_max = _solve_coupled_binary(wm, id_1, id_2, 1.0, _MOI.MAX_SENSE)

            if match_0_min == match_0_max
                append!(binary_variable_mapping, [((0.0, id_1), (match_0_min, id_2))])
            elseif match_1_min == match_1_max
                append!(binary_variable_mapping, [((1.0, id_1), (match_1_min, id_2))])
            end
        end
    end

    _revert_reduced_data!(data)
    return binary_variable_mapping
end


function _add_binary_cuts!(wm::AbstractWaterModel, binary_mapping)
    for nw in sort(collect(nw_ids(wm)))
        for mapping in binary_mapping
            index_1, index_2 = mapping[1][2], mapping[2][2]
            z_1 = var(wm, nw, index_1[3], index_1[4])
            z_2 = var(wm, nw, index_2[3], index_2[4])

            if mapping[1][1] == 0.0 && mapping[2][1] == 1.0
                JuMP.@constraint(wm.model, z_2 >= 1.0 - z_1)
            elseif mapping[1][1] == 1.0 && mapping[2][1] == 0.0
                JuMP.@constraint(wm.model, z_1 >= 1.0 - z_2)
            elseif mapping[1][1] == 1.0 && mapping[2][1] == 1.0
                JuMP.@constraint(wm.model, z_2 >= z_1)
            elseif mapping[1][1] == 0.0 && mapping[2][1] == 0.0
                JuMP.@constraint(wm.model, z_1 >= z_2)
            end
        end
    end
end


function _solve_coupled_binary(wm::AbstractWaterModel, index_1::Tuple, index_2::Tuple, fix_value::Float64, sense::_MOI.OptimizationSense)
    # Get the variable reference from the index tuple.
    z_1 = var(wm, index_1[1], index_1[3], index_1[4])
    z_2 = var(wm, index_2[1], index_2[3], index_2[4])

    # Fix the first variable.
    JuMP.fix(z_1, fix_value)

    # Optimize the variable (or affine expression) being tightened.
    JuMP.@objective(wm.model, sense, z_2)
    JuMP.optimize!(wm.model)
    termination_status = JuMP.termination_status(wm.model)

    if termination_status in [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL]
        candidate = JuMP.objective_value(wm.model) >= 0.5 ? 1.0 : 0.0
    else
        candidate = sense === _MOI.MIN_SENSE ? 0.0 : 1.0
    end

    # Unfix the first variable.
    JuMP.unfix(z_1)

    return candidate
end
