function correct_short_pipes!(data::Dict{String, <:Any})
    apply_wm!(_correct_short_pipes!, data; apply_to_subnetworks = true)
end


function correct_ne_short_pipes!(data::Dict{String, <:Any})
    apply_wm!(_correct_ne_short_pipes!, data; apply_to_subnetworks = true)
end


function _correct_short_pipes!(data::Dict{String, <:Any})
    capacity = calc_capacity_max(data)

    for (idx, short_pipe) in data["short_pipe"]
        _correct_status!(short_pipe)
        _correct_flow_direction!(short_pipe)
        _correct_short_pipe_flow_bounds!(short_pipe, capacity)
    end
end


function _correct_ne_short_pipes!(data::Dict{String, <:Any})
    capacity = calc_capacity_max(data)

    for (idx, ne_short_pipe) in data["ne_short_pipe"]
        _correct_status!(ne_short_pipe)
        _correct_flow_direction!(ne_short_pipe)
        _correct_short_pipe_flow_bounds!(ne_short_pipe, capacity)
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
    flow_min_dir = short_pipe["flow_direction"] == FLOW_DIRECTION_POSITIVE ? 0.0 : -Inf
    return max(-capacity, flow_min_dir, get(short_pipe, "flow_min", -Inf))
end


function _calc_short_pipe_flow_max(short_pipe::Dict{String, <:Any}, capacity::Float64)
    flow_max_dir = short_pipe["flow_direction"] == FLOW_DIRECTION_NEGATIVE ? 0.0 : Inf
    return min(capacity, flow_max_dir, get(short_pipe, "flow_max", Inf))
end


function set_short_pipe_warm_start!(data::Dict{String, <:Any})
    apply_wm!(_set_short_pipe_warm_start!, data)
end


function _set_short_pipe_warm_start!(data::Dict{String, <:Any})
    for short_pipe in values(data["short_pipe"])
        flow_mid = 0.5 * (short_pipe["flow_min"] + short_pipe["flow_max"])

        short_pipe["q_start"] = get(short_pipe, "q", flow_mid)
        short_pipe["qp_start"] = max(0.0, get(short_pipe, "q", flow_mid))
        short_pipe["qn_start"] = max(0.0, -get(short_pipe, "q", flow_mid))
    end
end


function set_ne_short_pipe_warm_start!(data::Dict{String, <:Any})
    apply_wm!(_set_ne_short_pipe_warm_start!, data)
end


function _set_ne_short_pipe_warm_start!(data::Dict{String, <:Any})
    for ne_short_pipe in values(data["ne_short_pipe"])
        flow_mid = 0.5 * (ne_short_pipe["flow_min"] + ne_short_pipe["flow_max"])

        ne_short_pipe["q_start"] = get(ne_short_pipe, "q", flow_mid)
        ne_short_pipe["qp_start"] = max(0.0, get(ne_short_pipe, "q", flow_mid))
        ne_short_pipe["qn_start"] = max(0.0, -get(ne_short_pipe, "q", flow_mid))
        ne_short_pipe["z_ne_short_pipe_start"] = abs(get(ne_short_pipe, "q", 0.0)) > 0.0 ? 1.0 : 0.0
    end
end
