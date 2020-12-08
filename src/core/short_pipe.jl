function correct_short_pipes!(data::Dict{String, <:Any})
    capacity = _calc_capacity_max(data)

    for (idx, short_pipe) in data["short_pipe"]
        _correct_short_pipe_flow_bounds!(short_pipe, capacity)
    end
end


function _correct_short_pipe_flow_bounds!(short_pipe::Dict{String, <:Any}, capacity::Float64)
    flow_min = _calc_short_pipe_flow_min(short_pipe, capacity)
    flow_max = _calc_short_pipe_flow_max(short_pipe, capacity)
    short_pipe["flow_min"], short_pipe["flow_max"] = flow_min, flow_max
    short_pipe["flow_min_forward"] = max(flow_min, get(short_pipe, "flow_min_forward", 0.0))
    short_pipe["flow_max_reverse"] = min(flow_max, get(short_pipe, "flow_max_reverse", 0.0))
end


function _calc_short_pipe_flow_min(short_pipe::Dict{String, <:Any}, capacity::Float64)
    flow_min_dir = short_pipe["flow_direction"] == POSITIVE ? 0.0 : -Inf
    return max(-capacity, flow_min_dir, get(short_pipe, "flow_min", -Inf))
end


function _calc_short_pipe_flow_max(short_pipe::Dict{String, <:Any}, capacity::Float64)
    flow_max_dir = short_pipe["flow_direction"] == NEGATIVE ? 0.0 : Inf
    return min(capacity, flow_max_dir, get(short_pipe, "flow_max", Inf))
end
