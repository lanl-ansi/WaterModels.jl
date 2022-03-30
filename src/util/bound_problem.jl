mutable struct BoundProblem
    sense::JuMP.OptimizationSense # Minimization or maximization.
    variable_to_tighten::_VariableIndex # Index of variable to tighten.
    variables_fix_one::Array{_VariableIndex} # Fix to one.
    variables_fix_zero::Array{_VariableIndex} # Fix to zero.
    bound_name::String # e.g., flow_min, flow_max, flow_min_forward
    bound::Float64 # The current bound of the variable.
    precision::Float64 # The desired precision of the variable bound.
    changed::Bool # Indicator specifying if the bound has changed.
    infeasible::Bool # Indicator specifying if the problem is infeasible.
end


function _var_exists(wm::AbstractWaterModel, nw_id::Int, var_symbol::Symbol, index::Int)
    if haskey(var(wm, nw_id), var_symbol)
        var_indices = var(wm, nw_id, var_symbol).axes[1]
        return index in var_indices
    else
        return false
    end
end


function _get_bound_problems_nodes(wm::AbstractWaterModel; kwargs...)
    return vcat(_get_bound_problems_nodes.(Ref(wm), nw_ids(wm); kwargs...)...)
end


function _get_bound_problems_nodes(wm::AbstractWaterModel, nw::Int; kwargs...)
    return vcat(_get_bound_problems_node.(Ref(wm), ids(wm, nw, :node), nw; kwargs...)...)
end


function _get_bound_problems_node(
    wm::AbstractWaterModel, i::Int, nw::Int)::Vector{BoundProblem}
    # Initialize the bound problem vector.
    bound_problems = Vector{BoundProblem}([])

    if _var_exists(wm, nw, :h, i)
        # Get the head variable at the node index.
        h_vid = _VariableIndex(nw, :node, :h, i)

        # Get the WaterModels portion of the data dictionary.
        wm_data = get_wm_data(wm.data)

        # Define the desired precision of the head bound.
        head_transform = _calc_head_per_unit_transform(wm_data)
        head_precision = head_transform(1.0e-3) # 1.0 mm.

        # Get the node corresponding to the index.
        node = ref(wm, nw, :node, i)

        # Get minimum head bounds from data and the variable and take the maximum.
        h_min_data = get(node, "h_min", -Inf)
        h_min_var = _get_lower_bound_from_index(wm, h_vid)
        h_min = max(h_min_data, h_min_var)

        # Get maximum head bounds from data and the variable and take the maximum.
        h_max_data = get(node, "h_max", Inf)
        h_max_var = _get_upper_bound_from_index(wm, h_vid)
        h_max = min(h_max_data, h_max_var)

        # Ensure the minimum and maximum bounds are sensible.
        @assert h_min <= h_max

        if h_min < h_max
            # Push a BoundProblem instance to minimize head.
            push!(bound_problems, BoundProblem(JuMP.MIN_SENSE, h_vid,
                [], [], "head_min", h_min, head_precision, true, false))

            # Push a BoundProblem instance to maximize head.
            push!(bound_problems, BoundProblem(JuMP.MAX_SENSE, h_vid,
                [], [], "head_max", h_max, head_precision, true, false))
        end
    end

    return bound_problems
end


function _get_bound_problems_tanks(wm::AbstractWaterModel; limit::Bool = false)
    return vcat(_get_bound_problems_tanks.(Ref(wm), nw_ids(wm); limit = limit)...)
end


function _get_bound_problems_tanks(wm::AbstractWaterModel, nw::Int; limit::Bool = false)
    return vcat(_get_bound_problems_tank.(Ref(wm), ids(wm, nw, :tank), nw; limit = limit)...)
end


function _get_bound_problems_tank(wm::AbstractWaterModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :q_tank) && i in [x for x in var(wm, nw, :q_tank).axes[1]]
        tank = ref(wm, nw, :tank, i)
        q_tank_vid = _VariableIndex(nw, :tank, :q_tank, i)
       
        wm_data = get_wm_data(wm.data)
        flow_transform = _calc_flow_per_unit_transform(wm_data)
        flow_precision = flow_transform(1.0e-5)

        q_tank_min = max(get(tank, "flow_min", -Inf), _get_lower_bound_from_index(wm, q_tank_vid))
        bp_min = BoundProblem(JuMP.MIN_SENSE, q_tank_vid, [], [], "flow_min", q_tank_min, flow_precision, true, false)

        q_tank_max = min(get(tank, "flow_max", Inf), _get_upper_bound_from_index(wm, q_tank_vid))
        bp_max = BoundProblem(JuMP.MAX_SENSE, q_tank_vid, [], [], "flow_max", q_tank_max, flow_precision, true, false)

        return Vector{BoundProblem}([bp_min, bp_max])
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems_pipes(wm::AbstractWaterModel; limit::Bool = false)
    return vcat(_get_bound_problems_pipes.(Ref(wm), nw_ids(wm); limit = limit)...)
end


function _get_bound_problems_pipes(wm::AbstractWaterModel, nw::Int; limit::Bool = false)
    return vcat(_get_bound_problems_pipe.(Ref(wm), ids(wm, nw, :pipe), nw; limit = limit)...)
end


function _get_bound_problems_pipe(wm::AbstractNCModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :q_pipe) && i in [x for x in var(wm, nw, :q_pipe).axes[1]]
        pipe = ref(wm, nw, :pipe, i)
        q_vid = _VariableIndex(nw, :pipe, :q_pipe, i)

        wm_data = get_wm_data(wm.data)
        flow_transform = _calc_flow_per_unit_transform(wm_data)
        flow_precision = flow_transform(1.0e-5)

        flow_min = max(get(pipe, "flow_min", -Inf), _get_lower_bound_from_index(wm, q_vid))
        bp_min = BoundProblem(JuMP.MIN_SENSE, q_vid, [], [], "flow_min", flow_min, flow_precision, true, false)

        flow_max = min(get(pipe, "flow_max", Inf), _get_upper_bound_from_index(wm, q_vid))
        bp_max = BoundProblem(JuMP.MAX_SENSE, q_vid, [], [], "flow_max", flow_max, flow_precision, true, false)

        return Vector{BoundProblem}([bp_min, bp_max])
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems_pipe(wm::AbstractNCDModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :q_pipe) && i in [x for x in var(wm, nw, :q_pipe).axes[1]]
        pipe = ref(wm, nw, :pipe, i)

        q_vid = _VariableIndex(nw, :pipe, :q_pipe, i)
        y_vid = _VariableIndex(nw, :pipe, :y_pipe, i)
        
        wm_data = get_wm_data(wm.data)
        flow_transform = _calc_flow_per_unit_transform(wm_data)
        flow_precision = flow_transform(1.0e-5)

        y_min = get(pipe, "y_min", 0.0)
        bp_y_min = BoundProblem(JuMP.MIN_SENSE, y_vid, [], [], "y_min", y_min, 1.0e-2, true, false)
       
        y_max = get(pipe, "y_max", 1.0)
        bp_y_max = BoundProblem(JuMP.MAX_SENSE, y_vid, [], [], "y_max", y_max, 1.0e-2, true, false)

        flow_min = max(get(pipe, "flow_min", -Inf), _get_lower_bound_from_index(wm, q_vid))
        flow_min = y_min == 1.0 ? max(0.0, flow_min) : flow_min
        bp_q_min = BoundProblem(JuMP.MIN_SENSE, q_vid, [], [], "flow_min", flow_min, flow_precision, true, false)
 
        flow_min_forward = get(pipe, "flow_min_forward", 0.0)
        flow_min_forward = y_max == 0.0 ? 0.0 : flow_min_forward
        bp_q_min_forward = BoundProblem(JuMP.MIN_SENSE, q_vid, [y_vid], [], "flow_min_forward", flow_min_forward, flow_precision, true, false)

        flow_max = min(get(pipe, "flow_max", Inf), _get_upper_bound_from_index(wm, q_vid))
        flow_max = y_max == 0.0 ? min(0.0, flow_max) : flow_max
        bp_q_max = BoundProblem(JuMP.MAX_SENSE, q_vid, [], [], "flow_max", flow_max, flow_precision, true, false)

        flow_max_reverse = get(pipe, "flow_max_reverse", 0.0)
        flow_max_reverse = y_min == 1.0 ? 0.0 : flow_max_reverse
        bp_q_max_reverse = BoundProblem(JuMP.MAX_SENSE, q_vid, [], [y_vid], "flow_max_reverse", flow_max_reverse, flow_precision, true, false)

        if limit
            return Vector{BoundProblem}([bp_q_min, bp_q_max])
        else
            return Vector{BoundProblem}([bp_q_min, bp_q_min_forward, bp_q_max, bp_q_max_reverse, bp_y_min, bp_y_max])
        end
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems_des_pipes(wm::AbstractWaterModel; limit::Bool = false)
    return vcat(_get_bound_problems_des_pipes.(Ref(wm), nw_ids(wm); limit = limit)...)
end


function _get_bound_problems_des_pipes(wm::AbstractWaterModel, nw::Int; limit::Bool = false)
    return vcat(_get_bound_problems_des_pipe.(Ref(wm), ids(wm, nw, :des_pipe), nw; limit = limit)...)
end


function _get_bound_problems_des_pipe(wm::AbstractNCDModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :q_des_pipe) && i in [x for x in var(wm, nw, :q_des_pipe).axes[1]]
        des_pipe = ref(wm, nw, :des_pipe, i)

        q_vid = _VariableIndex(nw, :des_pipe, :q_des_pipe, i)
        y_vid = _VariableIndex(nw, :des_pipe, :y_des_pipe, i)
        z_vid = _VariableIndex(nw, :des_pipe, :z_des_pipe, i)

        wm_data = get_wm_data(wm.data)
        flow_transform = _calc_flow_per_unit_transform(wm_data)
        flow_precision = flow_transform(1.0e-5)

        y_min = get(des_pipe, "y_min", 0.0)
        bp_y_min = BoundProblem(JuMP.MIN_SENSE, y_vid, [], [], "y_min", y_min, 1.0e-2, true, false)

        y_max = get(des_pipe, "y_max", 1.0)
        bp_y_max = BoundProblem(JuMP.MAX_SENSE, y_vid, [], [], "y_max", y_max, 1.0e-2, true, false)

        z_min = get(des_pipe, "z_min", 0.0)
        bp_z_min = BoundProblem(JuMP.MIN_SENSE, z_vid, [], [], "z_min", z_min, 1.0e-2, true, false)

        z_max = get(des_pipe, "z_max", 1.0)
        bp_z_max = BoundProblem(JuMP.MAX_SENSE, z_vid, [], [], "z_max", z_max, 1.0e-2, true, false)

        flow_min = max(get(des_pipe, "flow_min", -Inf), _get_lower_bound_from_index(wm, q_vid))
        flow_min = y_min == 1.0 ? max(0.0, flow_min) : flow_min
        flow_min = z_max == 0.0 ? 0.0 : flow_min
        bp_q_min = BoundProblem(JuMP.MIN_SENSE, q_vid, [], [], "flow_min", flow_min, flow_precision, true, false)

        flow_min_forward = get(des_pipe, "flow_min_forward", 0.0)
        flow_min_forward = y_max == 0.0 ? 0.0 : flow_min_forward
        flow_min_forward = z_max == 0.0 ? 0.0 : flow_min_forward
        bp_q_min_forward = BoundProblem(JuMP.MIN_SENSE, q_vid, [y_vid, z_vid], [], "flow_min_forward", flow_min_forward, flow_precision, true, false)

        flow_max = min(get(des_pipe, "flow_max", Inf), _get_upper_bound_from_index(wm, q_vid))
        flow_max = y_max == 0.0 ? min(0.0, flow_max) : flow_max
        flow_max = z_max == 0.0 ? 0.0 : flow_max
        bp_q_max = BoundProblem(JuMP.MAX_SENSE, q_vid, [], [], "flow_max", flow_max, flow_precision, true, false)

        flow_max_reverse = get(des_pipe, "flow_max_reverse", 0.0)
        flow_max_reverse = y_min == 1.0 ? 0.0 : flow_max_reverse
        flow_max_reverse = z_max == 0.0 ? 0.0 : flow_max_reverse
        bp_q_max_reverse = BoundProblem(JuMP.MAX_SENSE, q_vid, [z_vid], [y_vid], "flow_max_reverse", flow_max_reverse, flow_precision, true, false)

        if limit
            return Vector{BoundProblem}([bp_q_min_forward, bp_q_max_reverse])
        else
            return Vector{BoundProblem}([bp_q_min, bp_q_min_forward, bp_q_max, bp_q_max_reverse, bp_y_min, bp_y_max, bp_z_min, bp_z_max])
        end
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems_pumps(wm::AbstractWaterModel; limit::Bool = false)
    return vcat(_get_bound_problems_pumps.(Ref(wm), nw_ids(wm); limit = limit)...)
end


function _get_bound_problems_pumps(wm::AbstractWaterModel, nw::Int; limit::Bool = false)
    return vcat(_get_bound_problems_pump.(Ref(wm), ids(wm, nw, :pump), nw; limit = limit)...)
end


function _get_bound_problems_pump(wm::AbstractWaterModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :q_pump) && i in [x for x in var(wm, nw, :q_pump).axes[1]]
        pump = ref(wm, nw, :pump, i)

        q_vid = _VariableIndex(nw, :pump, :q_pump, i)
        z_vid = _VariableIndex(nw, :pump, :z_pump, i)

        wm_data = get_wm_data(wm.data)
        flow_transform = _calc_flow_per_unit_transform(wm_data)
        flow_precision = flow_transform(1.0e-5)

        z_min = get(pump, "z_min", 0.0)
        bp_z_min = BoundProblem(JuMP.MIN_SENSE, z_vid, [], [], "z_min", z_min, 1.0e-2, true, false)

        z_max = get(pump, "z_max", 1.0)
        bp_z_max = BoundProblem(JuMP.MAX_SENSE, z_vid, [], [], "z_max", z_max, 1.0e-2, true, false)

        flow_min_forward = get(pump, "flow_min_forward", 0.0)
        flow_min_forward = z_max == 0.0 ? 0.0 : flow_min_forward
        bp_min = BoundProblem(JuMP.MIN_SENSE, q_vid, [z_vid], [], "flow_min_forward", flow_min_forward, flow_precision, true, false)

        flow_max = min(get(pump, "flow_max", Inf), _get_upper_bound_from_index(wm, q_vid))
        flow_max = z_max == 0.0 ? 0.0 : flow_max
        bp_max = BoundProblem(JuMP.MAX_SENSE, q_vid, [], [], "flow_max", flow_max, flow_precision, true, false)

        return Vector{BoundProblem}([bp_min, bp_max, bp_z_min, bp_z_max])
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems_regulators(wm::AbstractWaterModel; limit::Bool = false)
    return vcat(_get_bound_problems_regulators.(Ref(wm), nw_ids(wm); limit = limit)...)
end


function _get_bound_problems_regulators(wm::AbstractWaterModel, nw::Int; limit::Bool = false)
    return vcat(_get_bound_problems_regulator.(Ref(wm), ids(wm, nw, :regulator), nw; limit = limit)...)
end


function _get_bound_problems_regulator(wm::AbstractWaterModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :q_regulator) && i in [x for x in var(wm, nw, :q_regulator).axes[1]]
        regulator = ref(wm, nw, :regulator, i)

        q_vid = _VariableIndex(nw, :regulator, :q_regulator, i)
        z_vid = _VariableIndex(nw, :regulator, :z_regulator, i)

        wm_data = get_wm_data(wm.data)
        flow_transform = _calc_flow_per_unit_transform(wm_data)
        flow_precision = flow_transform(1.0e-5)

        z_min = get(regulator, "z_min", 0.0)
        bp_z_min = BoundProblem(JuMP.MIN_SENSE, z_vid, [], [], "z_min", z_min, 1.0e-2, true, false)

        z_max = get(regulator, "z_max", 1.0)
        bp_z_max = BoundProblem(JuMP.MAX_SENSE, z_vid, [], [], "z_max", z_max, 1.0e-2, true, false)

        flow_min_forward = get(regulator, "flow_min_forward", 0.0)
        flow_min_forward = z_max == 0.0 ? 0.0 : flow_min_forward
        bp_min = BoundProblem(JuMP.MIN_SENSE, q_vid, [z_vid], [], "flow_min_forward", flow_min_forward, flow_precision, true, false)

        flow_max = min(get(regulator, "flow_max", Inf), _get_upper_bound_from_index(wm, q_vid))
        flow_max = z_max == 0.0 ? 0.0 : flow_max
        bp_max = BoundProblem(JuMP.MAX_SENSE, q_vid, [], [], "flow_max", flow_max, flow_precision, true, false)

        return Vector{BoundProblem}([bp_min, bp_max, bp_z_min, bp_z_max])
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems_short_pipes(wm::AbstractWaterModel; limit::Bool = false)
    return vcat(_get_bound_problems_short_pipes.(Ref(wm), nw_ids(wm); limit = limit)...)
end


function _get_bound_problems_short_pipes(wm::AbstractWaterModel, nw::Int; limit::Bool = false)
    return vcat(_get_bound_problems_short_pipe.(Ref(wm), ids(wm, nw, :short_pipe), nw; limit = limit)...)
end


function _get_bound_problems_short_pipe(wm::AbstractNCModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :q_short_pipe) && i in [x for x in var(wm, nw, :q_short_pipe).axes[1]]
        short_pipe = ref(wm, nw, :short_pipe, i)
        q_vid = _VariableIndex(nw, :short_pipe, :q_short_pipe, i)

        wm_data = get_wm_data(wm.data)
        flow_transform = _calc_flow_per_unit_transform(wm_data)
        flow_precision = flow_transform(1.0e-5)

        flow_min = max(get(short_pipe, "flow_min", -Inf), _get_lower_bound_from_index(wm, q_vid))
        bp_min = BoundProblem(JuMP.MIN_SENSE, q_vid, [], [], "flow_min", flow_min, flow_precision, true, false)

        flow_max = min(get(short_pipe, "flow_max", Inf), _get_lower_bound_from_index(wm, q_vid))
        bp_max = BoundProblem(JuMP.MAX_SENSE, q_vid, [], [], "flow_max", flow_max, flow_precision, true, false)

        return Vector{BoundProblem}([bp_min, bp_max])
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems_short_pipe(wm::AbstractNCDModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :q_short_pipe) && i in [x for x in var(wm, nw, :q_short_pipe).axes[1]]
        short_pipe = ref(wm, nw, :short_pipe, i)
        
        q_vid = _VariableIndex(nw, :short_pipe, :q_short_pipe, i)
        y_vid = _VariableIndex(nw, :short_pipe, :y_short_pipe, i)

        wm_data = get_wm_data(wm.data)
        flow_transform = _calc_flow_per_unit_transform(wm_data)
        flow_precision = flow_transform(1.0e-5)

        y_min = get(short_pipe, "y_min", 0.0)
        bp_y_min = BoundProblem(JuMP.MIN_SENSE, y_vid, [], [], "y_min", y_min, 1.0e-2, true, false)
       
        y_max = get(short_pipe, "y_max", 1.0)
        bp_y_max = BoundProblem(JuMP.MAX_SENSE, y_vid, [], [], "y_max", y_max, 1.0e-2, true, false)

        flow_min = max(get(short_pipe, "flow_min", -Inf), _get_lower_bound_from_index(wm, q_vid))
        flow_min = y_min == 1.0 ? max(0.0, flow_min) : flow_min
        bp_q_min = BoundProblem(JuMP.MIN_SENSE, q_vid, [], [], "flow_min", flow_min, flow_precision, true, false)
        
        flow_min_forward = get(short_pipe, "flow_min_forward", 0.0)
        flow_min_forward = y_max == 0.0 ? 0.0 : flow_min_forward
        bp_q_min_forward = BoundProblem(JuMP.MIN_SENSE, q_vid, [y_vid], [], "flow_min_forward", flow_min_forward, flow_precision, true, false)

        flow_max = min(get(short_pipe, "flow_max", Inf), _get_upper_bound_from_index(wm, q_vid))
        flow_max = y_max == 0.0 ? min(0.0, flow_max) : flow_max
        bp_q_max = BoundProblem(JuMP.MAX_SENSE, q_vid, [], [], "flow_max", flow_max, flow_precision, true, false)

        flow_max_reverse = get(short_pipe, "flow_max_reverse", 0.0)
        flow_max_reverse = y_min == 1.0 ? 0.0 : flow_max_reverse
        bp_q_max_reverse = BoundProblem(JuMP.MAX_SENSE, q_vid, [], [y_vid], "flow_max_reverse", flow_max_reverse, flow_precision, true, false)

        if limit
            return Vector{BoundProblem}([bp_q_min, bp_q_max])
        else
            return Vector{BoundProblem}([bp_q_min, bp_q_min_forward, bp_q_max, bp_q_max_reverse, bp_y_min, bp_y_max])
        end
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems_ne_short_pipes(wm::AbstractWaterModel; limit::Bool = false)
    return vcat(_get_bound_problems_ne_short_pipes.(Ref(wm), nw_ids(wm); limit = limit)...)
end


function _get_bound_problems_ne_short_pipes(wm::AbstractWaterModel, nw::Int; limit::Bool = false)
    return vcat(_get_bound_problems_ne_short_pipe.(Ref(wm), ids(wm, nw, :ne_short_pipe), nw; limit = limit)...)
end


function _get_bound_problems_ne_short_pipe(wm::AbstractNCModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :q_ne_short_pipe) && i in [x for x in var(wm, nw, :q_ne_short_pipe).axes[1]]
        ne_short_pipe = ref(wm, nw, :ne_short_pipe, i)

        q_vid = _VariableIndex(nw, :ne_short_pipe, :q_ne_short_pipe, i)
        z_vid = _VariableIndex(nw, :ne_short_pipe, :z_ne_short_pipe, i)

        wm_data = get_wm_data(wm.data)
        flow_transform = _calc_flow_per_unit_transform(wm_data)
        flow_precision = flow_transform(1.0e-5)

        z_min = get(ne_short_pipe, "z_min", 0.0)
        bp_z_min = BoundProblem(JuMP.MIN_SENSE, z_vid, [], [], "z_min", z_min, 1.0e-2, true, false)

        z_max = get(ne_short_pipe, "z_max", 1.0)
        bp_z_max = BoundProblem(JuMP.MAX_SENSE, z_vid, [], [], "z_max", z_max, 1.0e-2, true, false)

        flow_min = max(get(ne_short_pipe, "flow_min", -Inf), _get_lower_bound_from_index(wm, q_vid))
        flow_min = z_max == 0.0 ? 0.0 : flow_min
        bp_min = BoundProblem(JuMP.MIN_SENSE, q_vid, [], [], "flow_min", flow_min, flow_precision, true, false)

        flow_max = min(get(ne_short_pipe, "flow_max", Inf), _get_upper_bound_from_index(wm, q_vid))
        flow_max = z_max == 0.0 ? 0.0 : flow_max
        bp_max = BoundProblem(JuMP.MAX_SENSE, q_vid, [], [], "flow_max", flow_max, flow_precision, true, false)

        return Vector{BoundProblem}([bp_min, bp_max, bp_z_min, bp_z_max])
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems_ne_short_pipe(wm::AbstractNCDModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :q_ne_short_pipe) && i in [x for x in var(wm, nw, :q_ne_short_pipe).axes[1]]
        ne_short_pipe = ref(wm, nw, :ne_short_pipe, i)

        q_vid = _VariableIndex(nw, :ne_short_pipe, :q_ne_short_pipe, i)
        y_vid = _VariableIndex(nw, :ne_short_pipe, :y_ne_short_pipe, i)
        z_vid = _VariableIndex(nw, :ne_short_pipe, :z_ne_short_pipe, i)

        wm_data = get_wm_data(wm.data)
        flow_transform = _calc_flow_per_unit_transform(wm_data)
        flow_precision = flow_transform(1.0e-5)

        y_min = get(ne_short_pipe, "y_min", 0.0)
        bp_y_min = BoundProblem(JuMP.MIN_SENSE, y_vid, [], [], "y_min", y_min, 1.0e-2, true, false)

        y_max = get(ne_short_pipe, "y_max", 1.0)
        bp_y_max = BoundProblem(JuMP.MAX_SENSE, y_vid, [], [], "y_max", y_max, 1.0e-2, true, false)

        z_min = get(ne_short_pipe, "z_min", 0.0)
        bp_z_min = BoundProblem(JuMP.MIN_SENSE, z_vid, [], [], "z_min", z_min, 1.0e-2, true, false)

        z_max = get(ne_short_pipe, "z_max", 1.0)
        bp_z_max = BoundProblem(JuMP.MAX_SENSE, z_vid, [], [], "z_max", z_max, 1.0e-2, true, false)

        flow_min = max(get(ne_short_pipe, "flow_min", -Inf), _get_lower_bound_from_index(wm, q_vid))
        flow_min = y_min == 1.0 ? max(0.0, flow_min) : flow_min
        flow_min = z_max == 0.0 ? 0.0 : flow_min
        bp_q_min = BoundProblem(JuMP.MIN_SENSE, q_vid, [], [], "flow_min", flow_min, flow_precision, true, false)
        
        flow_min_forward = get(ne_short_pipe, "flow_min_forward", 0.0)
        flow_min_forward = y_max == 0.0 ? 0.0 : flow_min_forward
        flow_min_forward = z_max == 0.0 ? 0.0 : flow_min_forward
        bp_q_min_forward = BoundProblem(JuMP.MIN_SENSE, q_vid, [y_vid, z_vid], [], "flow_min_forward", flow_min_forward, flow_precision, true, false)

        flow_max = min(get(ne_short_pipe, "flow_max", Inf), _get_upper_bound_from_index(wm, q_vid))
        flow_max = y_max == 0.0 ? min(0.0, flow_max) : flow_max
        flow_max = z_max == 0.0 ? 0.0 : flow_max
        bp_q_max = BoundProblem(JuMP.MAX_SENSE, q_vid, [], [], "flow_max", flow_max, flow_precision, true, false)

        flow_max_reverse = get(ne_short_pipe, "flow_max_reverse", 0.0)
        flow_max_reverse = y_min == 1.0 ? 0.0 : flow_max_reverse
        flow_max_reverse = z_max == 0.0 ? 0.0 : flow_max_reverse
        bp_q_max_reverse = BoundProblem(JuMP.MAX_SENSE, q_vid, [z_vid], [y_vid], "flow_max_reverse", flow_max_reverse, flow_precision, true, false)

        if limit
            return Vector{BoundProblem}([bp_q_min, bp_q_max])
        else
            return Vector{BoundProblem}([bp_q_min, bp_q_min_forward, bp_q_max, bp_q_max_reverse, bp_y_min, bp_y_max, bp_z_min, bp_z_max])
        end
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems_valves(wm::AbstractWaterModel; limit::Bool = false)
    return vcat(_get_bound_problems_valves.(Ref(wm), nw_ids(wm); limit = limit)...)
end


function _get_bound_problems_valves(wm::AbstractWaterModel, nw::Int; limit::Bool = false)
    return vcat(_get_bound_problems_valve.(Ref(wm), ids(wm, nw, :valve), nw; limit = limit)...)
end


function _get_bound_problems_valve(wm::AbstractNCModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :q_valve) && i in [x for x in var(wm, nw, :q_valve).axes[1]]
        valve = ref(wm, nw, :valve, i)

        q_vid = _VariableIndex(nw, :valve, :q_valve, i)
        z_vid = _VariableIndex(nw, :valve, :z_valve, i)

        wm_data = get_wm_data(wm.data)
        flow_transform = _calc_flow_per_unit_transform(wm_data)
        flow_precision = flow_transform(1.0e-5)

        z_min = get(valve, "z_min", 0.0)
        bp_z_min = BoundProblem(JuMP.MIN_SENSE, z_vid, [], [], "z_min", z_min, 1.0e-2, true, false)

        z_max = get(valve, "z_max", 1.0)
        bp_z_max = BoundProblem(JuMP.MAX_SENSE, z_vid, [], [], "z_max", z_max, 1.0e-2, true, false)

        flow_min = max(get(valve, "flow_min", -Inf), _get_lower_bound_from_index(wm, q_vid))
        flow_min = z_max == 0.0 ? 0.0 : flow_min
        bp_min = BoundProblem(JuMP.MIN_SENSE, q_vid, [], [], "flow_min", flow_min, flow_precision, true, false)

        flow_max = min(get(valve, "flow_max", Inf), _get_upper_bound_from_index(wm, q_vid))
        flow_max = z_max == 0.0 ? 0.0 : flow_max
        bp_max = BoundProblem(JuMP.MAX_SENSE, q_vid, [], [], "flow_max", flow_max, flow_precision, true, false)

        return Vector{BoundProblem}([bp_min, bp_max, bp_z_min, bp_z_max])
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems_valve(wm::AbstractNCDModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :q_valve) && i in [x for x in var(wm, nw, :q_valve).axes[1]]
        valve = ref(wm, nw, :valve, i)

        q_vid = _VariableIndex(nw, :valve, :q_valve, i)
        y_vid = _VariableIndex(nw, :valve, :y_valve, i)
        z_vid = _VariableIndex(nw, :valve, :z_valve, i)

        wm_data = get_wm_data(wm.data)
        flow_transform = _calc_flow_per_unit_transform(wm_data)
        flow_precision = flow_transform(1.0e-5)

        y_min = get(valve, "y_min", 0.0)
        bp_y_min = BoundProblem(JuMP.MIN_SENSE, y_vid, [], [], "y_min", y_min, 1.0e-2, true, false)

        y_max = get(valve, "y_max", 1.0)
        bp_y_max = BoundProblem(JuMP.MAX_SENSE, y_vid, [], [], "y_max", y_max, 1.0e-2, true, false)

        z_min = get(valve, "z_min", 0.0)
        bp_z_min = BoundProblem(JuMP.MIN_SENSE, z_vid, [], [], "z_min", z_min, 1.0e-2, true, false)

        z_max = get(valve, "z_max", 1.0)
        bp_z_max = BoundProblem(JuMP.MAX_SENSE, z_vid, [], [], "z_max", z_max, 1.0e-2, true, false)

        flow_min = max(get(valve, "flow_min", -Inf), _get_lower_bound_from_index(wm, q_vid))
        flow_min = y_min == 1.0 ? max(0.0, flow_min) : flow_min
        flow_min = z_max == 0.0 ? 0.0 : flow_min
        bp_q_min = BoundProblem(JuMP.MIN_SENSE, q_vid, [], [], "flow_min", flow_min, flow_precision, true, false)
        
        flow_min_forward = get(valve, "flow_min_forward", 0.0)
        flow_min_forward = y_max == 0.0 ? 0.0 : flow_min_forward
        flow_min_forward = z_max == 0.0 ? 0.0 : flow_min_forward
        bp_q_min_forward = BoundProblem(JuMP.MIN_SENSE, q_vid, [y_vid, z_vid], [], "flow_min_forward", flow_min_forward, flow_precision, true, false)

        flow_max = min(get(valve, "flow_max", Inf), _get_upper_bound_from_index(wm, q_vid))
        flow_max = y_max == 0.0 ? min(0.0, flow_max) : flow_max
        flow_max = z_max == 0.0 ? 0.0 : flow_max
        bp_q_max = BoundProblem(JuMP.MAX_SENSE, q_vid, [], [], "flow_max", flow_max, flow_precision, true, false)

        flow_max_reverse = get(valve, "flow_max_reverse", 0.0)
        flow_max_reverse = y_min == 1.0 ? 0.0 : flow_max_reverse
        flow_max_reverse = z_max == 0.0 ? 0.0 : flow_max_reverse
        bp_q_max_reverse = BoundProblem(JuMP.MAX_SENSE, q_vid, [z_vid], [y_vid], "flow_max_reverse", flow_max_reverse, flow_precision, true, false)

        if limit
            return Vector{BoundProblem}([bp_q_min, bp_q_max])
        else
            return Vector{BoundProblem}([bp_q_min, bp_q_min_forward, bp_q_max, bp_q_max_reverse, bp_y_min, bp_y_max, bp_z_min, bp_z_max])
        end
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems(wm::AbstractWaterModel; limit::Bool = false)::Vector{BoundProblem}
    # Create the sets of bound-tightening problems.
    bps_node = _get_bound_problems_nodes(wm)
    bps_tank = _get_bound_problems_tanks(wm; limit = limit)
    bps_pipe = _get_bound_problems_pipes(wm; limit = limit)
    bps_des_pipe = _get_bound_problems_des_pipes(wm; limit = limit)
    bps_pump = _get_bound_problems_pumps(wm; limit = limit)
    bps_regulator = _get_bound_problems_regulators(wm; limit = limit)
    bps_short_pipe = _get_bound_problems_short_pipes(wm; limit = limit)
    bps_ne_short_pipe = _get_bound_problems_ne_short_pipes(wm; limit = limit)
    bps_valve = _get_bound_problems_valves(wm; limit = limit)

    # Return the concatenation of all bound problem vectors.
    return vcat(bps_pump, bps_valve, bps_regulator, bps_des_pipe,
        bps_pipe, bps_short_pipe, bps_ne_short_pipe, bps_node, bps_tank)
end
