mutable struct BoundProblem
    sense::_MOI.OptimizationSense # Minimization or maximization.
    variable_to_tighten::_VariableIndex # Index of variable to tighten.
    variables_fix_one::Array{_VariableIndex} # Fix to one.
    variables_fix_zero::Array{_VariableIndex} # Fix to zero.
    bound_name::String # e.g., flow_min, flow_max, flow_min_forward
    bound::Float64 # The current bound of the variable.
    precision::Float64 # The precision of the variable bound.
    changed::Bool # Indicator specifying if the bound has changed.
end


function _get_bound_problems_nodes(wm::AbstractWaterModel; limit::Bool = false)
    return vcat(_get_bound_problems_nodes.(Ref(wm), nw_ids(wm); limit = limit)...)
end


function _get_bound_problems_nodes(wm::AbstractWaterModel, nw::Int; limit::Bool = false)
    return vcat(_get_bound_problems_node.(Ref(wm), ids(wm, nw, :node), nw; limit = limit)...)
end


function _get_bound_problems_node(wm::AbstractWaterModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :h) && i in [x[1] for x in eachindex(var(wm, nw, :h))]
        h_vid = _VariableIndex(nw, :node, :h, i)

        h_min = _get_lower_bound_from_index(wm, h_vid)
        bp_min = BoundProblem(_MOI.MIN_SENSE, h_vid, [],
            [], "head_min", h_min, 1.0e-2, true)

        h_max = _get_upper_bound_from_index(wm, h_vid)
        bp_max = BoundProblem(_MOI.MAX_SENSE, h_vid, [],
            [], "head_max", h_max, 1.0e-2, true)

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
    if haskey(var(wm, nw), :q_pipe) && i in [x[1] for x in eachindex(var(wm, nw, :q_pipe))]
        q_vid = _VariableIndex(nw, :pipe, :q_pipe, i)

        flow_min = _get_lower_bound_from_index(wm, q_vid)
        bp_min = BoundProblem(_MOI.MIN_SENSE, q_vid, [],
            [], "flow_min", flow_min, 1.0e-3, true)

        flow_max = _get_lower_bound_from_index(wm, q_vid)
        bp_max = BoundProblem(_MOI.MAX_SENSE, q_vid, [],
            [], "flow_max", flow_max, 1.0e-3, true)

        return Vector{BoundProblem}([bp_min, bp_max])
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems_pipe(wm::AbstractNCDModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :q_pipe) && i in [x[1] for x in eachindex(var(wm, nw, :q_pipe))]
        q_vid = _VariableIndex(nw, :pipe, :q_pipe, i)
        y_vid = _VariableIndex(nw, :pipe, :y_pipe, i)

        flow_min = _get_lower_bound_from_index(wm, q_vid)
        bp_q_min = BoundProblem(_MOI.MIN_SENSE, q_vid, [],
            [], "flow_min", flow_min, 1.0e-3, true)
        
        flow_min_forward = get(ref(wm, q_vid.network_index,
            :pipe)[i], "flow_min_forward", 0.0)
        bp_q_min_forward = BoundProblem(_MOI.MIN_SENSE, q_vid, [y_vid],
            [], "flow_min_forward", flow_min_forward, 1.0e-3, true)

        flow_max = _get_upper_bound_from_index(wm, q_vid)
        bp_q_max = BoundProblem(_MOI.MAX_SENSE, q_vid, [],
            [], "flow_max", flow_max, 1.0e-3, true)

        flow_max_reverse = get(ref(wm, q_vid.network_index,
            :pipe)[i], "flow_max_reverse", 0.0)
        bp_q_max_reverse = BoundProblem(_MOI.MAX_SENSE, q_vid, [],
            [y_vid], "flow_max_reverse", flow_max_reverse, 1.0e-3, true)

        bp_y_min = BoundProblem(_MOI.MIN_SENSE, y_vid, [],
            [], "y_min", 0.0, 1.0e-2, true)
        bp_y_max = BoundProblem(_MOI.MAX_SENSE, y_vid, [],
            [], "y_max", 1.0, 1.0e-2, true)

        if limit
            return Vector{BoundProblem}([bp_q_min, bp_q_max])
        else
            return Vector{BoundProblem}([bp_q_min, bp_q_min_forward,
                bp_q_max, bp_q_max_reverse, bp_y_min, bp_y_max])
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
    if haskey(var(wm, nw), :q_pump) && i in [x[1] for x in eachindex(var(wm, nw, :q_pump))]
        q_vid = _VariableIndex(nw, :pump, :q_pump, i)
        z_vid = _VariableIndex(nw, :pump, :z_pump, i)

        flow_min_forward = get(ref(wm, q_vid.network_index,
            :pump)[i], "flow_min_forward", 0.0)
        bp_min = BoundProblem(_MOI.MIN_SENSE, q_vid, [z_vid],
            [], "flow_min_forward", flow_min_forward, 1.0e-3, true)

        flow_max = _get_upper_bound_from_index(wm, q_vid)
        bp_max = BoundProblem(_MOI.MAX_SENSE, q_vid, [z_vid],
            [], "flow_max", flow_max, 1.0e-3, true)

        bp_z_min = BoundProblem(_MOI.MIN_SENSE, z_vid, [],
            [], "z_min", 0.0, 1.0e-2, true)
        bp_z_max = BoundProblem(_MOI.MAX_SENSE, z_vid, [],
            [], "z_max", 1.0, 1.0e-2, true)

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
    if haskey(var(wm, nw), :q_regulator) && i in [x[1] for x in eachindex(var(wm, nw, :q_regulator))]
        q_vid = _VariableIndex(nw, :regulator, :q_regulator, i)
        z_vid = _VariableIndex(nw, :regulator, :z_regulator, i)

        flow_min_forward = get(ref(wm, q_vid.network_index,
            :regulator)[i], "flow_min_forward", 0.0)
        bp_min = BoundProblem(_MOI.MIN_SENSE, q_vid, [z_vid],
            [], "flow_min_forward", flow_min_forward, 1.0e-3, true)

        flow_max = _get_upper_bound_from_index(wm, q_vid)
        bp_max = BoundProblem(_MOI.MAX_SENSE, q_vid, [z_vid],
            [], "flow_max", flow_max, 1.0e-3, true)

        bp_z_min = BoundProblem(_MOI.MIN_SENSE, z_vid, [],
            [], "z_min", 0.0, 1.0e-2, true)
        bp_z_max = BoundProblem(_MOI.MAX_SENSE, z_vid, [],
            [], "z_max", 1.0, 1.0e-2, true)

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
    if haskey(var(wm, nw), :q_short_pipe) && i in [x[1] for x in eachindex(var(wm, nw, :q_short_pipe))]
        q_vid = _VariableIndex(nw, :short_pipe, :q_short_pipe, i)

        flow_min = _get_lower_bound_from_index(wm, q_vid)
        bp_min = BoundProblem(_MOI.MIN_SENSE, q_vid, [],
            [], "flow_min", flow_min, 1.0e-3, true)

        flow_max = _get_lower_bound_from_index(wm, q_vid)
        bp_max = BoundProblem(_MOI.MAX_SENSE, q_vid, [],
            [], "flow_max", flow_max, 1.0e-3, true)

        return Vector{BoundProblem}([bp_min, bp_max])
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems_short_pipe(wm::AbstractNCDModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :q_short_pipe) && i in [x[1] for x in eachindex(var(wm, nw, :q_short_pipe))]
        q_vid = _VariableIndex(nw, :short_pipe, :q_short_pipe, i)
        y_vid = _VariableIndex(nw, :short_pipe, :y_short_pipe, i)

        flow_min = _get_lower_bound_from_index(wm, q_vid)
        bp_q_min = BoundProblem(_MOI.MIN_SENSE, q_vid, [],
            [], "flow_min", flow_min, 1.0e-3, true)
        
        flow_min_forward = get(ref(wm, q_vid.network_index,
            :short_pipe)[i], "flow_min_forward", 0.0)
        bp_q_min_forward = BoundProblem(_MOI.MIN_SENSE, q_vid, [y_vid],
            [], "flow_min_forward", flow_min_forward, 1.0e-3, true)

        flow_max = _get_upper_bound_from_index(wm, q_vid)
        bp_q_max = BoundProblem(_MOI.MAX_SENSE, q_vid, [],
            [], "flow_max", flow_max, 1.0e-3, true)

        flow_max_reverse = get(ref(wm, q_vid.network_index,
            :short_pipe)[i], "flow_max_reverse", 0.0)
        bp_q_max_reverse = BoundProblem(_MOI.MAX_SENSE, q_vid, [],
            [y_vid], "flow_max_reverse", flow_max_reverse, 1.0e-3, true)

        bp_y_min = BoundProblem(_MOI.MIN_SENSE, y_vid, [],
            [], "y_min", 0.0, 1.0e-2, true)
        bp_y_max = BoundProblem(_MOI.MAX_SENSE, y_vid, [],
            [], "y_max", 1.0, 1.0e-2, true)

        if limit
            return Vector{BoundProblem}([bp_q_min, bp_q_max])
        else
            return Vector{BoundProblem}([bp_q_min, bp_q_min_forward,
                bp_q_max, bp_q_max_reverse, bp_y_min, bp_y_max])
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
    if haskey(var(wm, nw), :q_valve) && i in [x[1] for x in eachindex(var(wm, nw, :q_valve))]
        q_vid = _VariableIndex(nw, :valve, :q_valve, i)
        z_vid = _VariableIndex(nw, :valve, :z_valve, i)

        flow_min = _get_lower_bound_from_index(wm, q_vid)
        bp_min = BoundProblem(_MOI.MIN_SENSE, q_vid, [],
            [], "flow_min", flow_min, 1.0e-3, true)

        flow_max = _get_lower_bound_from_index(wm, q_vid)
        bp_max = BoundProblem(_MOI.MAX_SENSE, q_vid, [],
            [], "flow_max", flow_max, 1.0e-3, true)

        bp_z_min = BoundProblem(_MOI.MIN_SENSE, z_vid, [],
            [], "z_min", 0.0, 1.0e-2, true)
        bp_z_max = BoundProblem(_MOI.MAX_SENSE, z_vid, [],
            [], "z_max", 1.0, 1.0e-2, true)

        return Vector{BoundProblem}([bp_min, bp_max, bp_z_min, bp_z_max])
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems_valve(wm::AbstractNCDModel, i::Int, nw::Int; limit::Bool = false)
    if haskey(var(wm, nw), :q_valve) && i in [x[1] for x in eachindex(var(wm, nw, :q_valve))]
        q_vid = _VariableIndex(nw, :valve, :q_valve, i)
        y_vid = _VariableIndex(nw, :valve, :y_valve, i)
        z_vid = _VariableIndex(nw, :valve, :z_valve, i)

        flow_min = _get_lower_bound_from_index(wm, q_vid)
        bp_q_min = BoundProblem(_MOI.MIN_SENSE, q_vid, [],
            [], "flow_min", flow_min, 1.0e-3, true)
        
        flow_min_forward = get(ref(wm, q_vid.network_index,
            :valve)[i], "flow_min_forward", 0.0)
        bp_q_min_forward = BoundProblem(_MOI.MIN_SENSE, q_vid, [y_vid, z_vid],
            [], "flow_min_forward", flow_min_forward, 1.0e-3, true)

        flow_max = _get_upper_bound_from_index(wm, q_vid)
        bp_q_max = BoundProblem(_MOI.MAX_SENSE, q_vid, [],
            [], "flow_max", flow_max, 1.0e-3, true)

        flow_max_reverse = get(ref(wm, q_vid.network_index,
            :valve)[i], "flow_max_reverse", 0.0)
        bp_q_max_reverse = BoundProblem(_MOI.MAX_SENSE, q_vid, [z_vid],
            [y_vid], "flow_max_reverse", flow_max_reverse, 1.0e-3, true)

        bp_y_min = BoundProblem(_MOI.MIN_SENSE, y_vid, [],
            [], "y_min", 0.0, 1.0e-2, true)
        bp_y_max = BoundProblem(_MOI.MAX_SENSE, y_vid, [],
            [], "y_max", 1.0, 1.0e-2, true)
        bp_z_min = BoundProblem(_MOI.MIN_SENSE, z_vid, [],
            [], "z_min", 0.0, 1.0e-2, true)
        bp_z_max = BoundProblem(_MOI.MAX_SENSE, z_vid, [],
            [], "z_max", 1.0, 1.0e-2, true)

        if limit
            return Vector{BoundProblem}([bp_q_min, bp_q_max])
        else
            return Vector{BoundProblem}([bp_q_min, bp_q_min_forward, bp_q_max,
                bp_q_max_reverse, bp_y_min, bp_y_max, bp_z_min, bp_z_max])
        end
    else
        return Vector{BoundProblem}([])
    end
end


function _get_bound_problems(wm::AbstractWaterModel; limit::Bool = false)::Vector{BoundProblem}
    # Create the sets of bound-tightening problems.
    bps_node = _get_bound_problems_nodes(wm; limit = limit)
    bps_pipe = _get_bound_problems_pipes(wm; limit = limit)
    bps_pump = _get_bound_problems_pumps(wm; limit = limit)
    bps_regulator = _get_bound_problems_regulators(wm; limit = limit)
    bps_short_pipe = _get_bound_problems_short_pipes(wm; limit = limit)
    bps_valve = _get_bound_problems_valves(wm; limit = limit)
    return vcat(bps_pump, bps_valve, bps_regulator,
        bps_pipe, bps_short_pipe, bps_node)
end