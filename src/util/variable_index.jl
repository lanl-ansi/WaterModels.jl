mutable struct _VariableIndex
    network_index::Int
    component_type::Symbol
    variable_symbol::Symbol
    component_index::Int
end


function _get_variable_from_index(wm::AbstractWaterModel, index::_VariableIndex)
    return var(wm, index.network_index, index.variable_symbol, index.component_index)
end


function _model_has_variable_index(wm::AbstractWaterModel, index::_VariableIndex)
    return index.component_index in ids(wm, index.network_index, index.component_type)
end


function _get_lower_bound_from_index(wm::AbstractWaterModel, index::_VariableIndex)
    v = _get_variable_from_index(wm, index)

    if typeof(v) === JuMP.VariableRef
        return JuMP.is_binary(v) ? 0.0 : JuMP.lower_bound(v)
    elseif typeof(v) == JuMP.AffExpr
        neg = collect(keys(filter(x -> x[2] < 0.0, v.terms)))
        pos = collect(keys(filter(x -> x[2] >= 0.0, v.terms)))
        return sum(JuMP.lower_bound.(pos)) - sum(JuMP.upper_bound.(neg))
    end
end


function _get_upper_bound_from_index(wm::AbstractWaterModel, index::_VariableIndex)
    v = _get_variable_from_index(wm, index)

    if typeof(v) === JuMP.VariableRef
        return JuMP.is_binary(v) ? 1.0 : JuMP.upper_bound(v)
    elseif typeof(v) == JuMP.AffExpr
        neg = collect(keys(filter(x -> x[2] < 0.0, v.terms)))
        pos = collect(keys(filter(x -> x[2] >= 0.0, v.terms)))
        return sum(JuMP.upper_bound.(pos)) - sum(JuMP.lower_bound.(neg))
    end
end


function _get_indicator_variable_indices(wm::AbstractWaterModel; nw::Int=nw_id_default)
    vars = Array{_VariableIndex, 1}()

    for comp_type in [:pump, :regulator, :valve]
        for comp_id in ids(wm, nw, comp_type)
            v_sym = Symbol("z_" * string(comp_type))
            append!(vars, [_VariableIndex(nw, comp_type, v_sym, comp_id)])
        end
    end

    return vars
end


function _get_flow_variable_indices(wm::AbstractWaterModel; nw::Int=nw_id_default)
    vars = Array{_VariableIndex, 1}()

    for comp_type in vcat(_LINK_COMPONENTS, "tank")
        for comp_id in ids(wm, nw, Symbol(comp_type))
            v_sym = Symbol("q_" * comp_type)
            append!(vars, [_VariableIndex(nw, Symbol(comp_type), v_sym, comp_id)])
        end
    end

    return vars
end


function _get_head_variable_indices(wm::AbstractWaterModel; nw::Int=nw_id_default)
    vars = Array{_VariableIndex, 1}()

    for i in ids(wm, nw, :node)
        append!(vars, [_VariableIndex(nw, :node, :h, i)])
    end

    return vars
end 


function _get_direction_variable_indices(wm::AbstractNCDModel; nw::Int=nw_id_default)
    vars = Array{_VariableIndex, 1}()

    for comp_type in _LINK_COMPONENTS
        for comp_id in ids(wm, nw, Symbol(comp_type))
            v_sym = Symbol("y_" * comp_type)
            append!(vars, [_VariableIndex(nw, Symbol(comp_type), v_sym, comp_id)])
        end
    end

    return vars
end


function _get_binary_variable_indices(wm::AbstractWaterModel; nw::Int=nw_id_default)
    return _get_indicator_variable_indices(wm; nw=nw)
end


function _get_binary_variable_indices(wm::AbstractNCDModel; nw::Int=nw_id_default)
    indicator_vars = _get_indicator_variable_indices(wm; nw=nw)
    direction_vars = _get_direction_variable_indices(wm; nw=nw)
    return vcat(indicator_vars, direction_vars)
end
