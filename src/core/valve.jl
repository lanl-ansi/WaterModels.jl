function correct_valves!(data::Dict{String, <:Any})
    apply_wm!(_correct_valves!, data; apply_to_subnetworks = true)
end


function _correct_valves!(data::Dict{String, <:Any})
    capacity = calc_capacity_max(data)

    for (idx, valve) in data["valve"]
        _correct_status!(valve)
        _correct_flow_direction!(valve)
        _correct_valve_flow_bounds!(valve, capacity)
    end
end


function _correct_valve_flow_bounds!(valve::Dict{String, <:Any}, capacity::Float64)
    flow_min = _calc_valve_flow_min(valve, capacity)
    flow_max = _calc_valve_flow_max(valve, capacity)
    valve["flow_min"], valve["flow_max"] = flow_min, flow_max
    valve["flow_min_forward"] = max(flow_min, get(valve, "flow_min_forward", 0.0))
    valve["flow_max_reverse"] = min(flow_max, get(valve, "flow_max_reverse", 0.0))

    if get(valve, "y_min", 0.0) == 1.0
        valve["flow_min"] = max(0.0, valve["flow_min"])
        valve["flow_min_forward"] = max(0.0, valve["flow_min_forward"])
        valve["flow_max_reverse"] = 0.0
        valve["flow_max"] = max(0.0, valve["flow_max"])
    elseif get(valve, "y_max", 1.0) == 0.0
        valve["flow_min"] = min(0.0, valve["flow_min"])
        valve["flow_min_forward"] = 0.0
        valve["flow_max_reverse"] = min(0.0, valve["flow_max_reverse"])
        valve["flow_max"] = min(0.0, valve["flow_max"])
    end

    if get(valve, "z_max", 1.0) == 0.0
        valve["flow_min"] = 0.0
        valve["flow_max"] = 0.0
        valve["flow_min_forward"] = 0.0
        valve["flow_max_reverse"] = 0.0
    end

    valve["flow_max"] = max(valve["flow_min"], valve["flow_max"])
    valve["flow_min"] = min(valve["flow_min"], valve["flow_max"])
    valve["flow_min_forward"] = max(valve["flow_min_forward"], valve["flow_min"])
    valve["flow_min_forward"] = min(valve["flow_min_forward"], valve["flow_max"])
    valve["flow_max_reverse"] = min(valve["flow_max_reverse"], valve["flow_max"])
    valve["flow_max_reverse"] = max(valve["flow_max_reverse"], valve["flow_min"])

    @assert valve["flow_min"] <= valve["flow_max"]
    @assert get(valve, "flow_min_forward", 0.0) <= max(0.0, valve["flow_max"])
    @assert min(0.0, valve["flow_min"]) <= get(valve, "flow_max_reverse", 0.0)
end


function _calc_valve_flow_min(valve::Dict{String, <:Any}, capacity::Float64)
    flow_min_dir = valve["flow_direction"] == FLOW_DIRECTION_POSITIVE ? 0.0 : -Inf
    return max(-capacity, flow_min_dir, get(valve, "flow_min", -Inf))
end


function _calc_valve_flow_max(valve::Dict{String, <:Any}, capacity::Float64)
    flow_max_dir = valve["flow_direction"] == FLOW_DIRECTION_NEGATIVE ? 0.0 : Inf
    return min(capacity, flow_max_dir, get(valve, "flow_max", Inf))
end


function _relax_valves!(data::Dict{String,<:Any})
    if !_IM.ismultinetwork(data)
        if haskey(data, "time_series") && haskey(data["time_series"], "valve")
            ts = data["time_series"]["valve"]
            valves = values(filter(x -> x.first in keys(ts), data["valve"]))
            map(x -> x["flow_min"] = minimum(ts[string(x["index"])]["flow_min"]), valves)
            map(x -> x["flow_min_forward"] = minimum(ts[string(x["index"])]["flow_min_forward"]), valves)
            map(x -> x["flow_max"] = maximum(ts[string(x["index"])]["flow_max"]), valves)
            map(x -> x["flow_max_reverse"] = maximum(ts[string(x["index"])]["flow_max_reverse"]), valves)
            map(x -> x["y_min"] = minimum(ts[string(x["index"])]["y_min"]), valves)
            map(x -> x["y_max"] = maximum(ts[string(x["index"])]["y_max"]), valves)
            map(x -> x["z_min"] = minimum(ts[string(x["index"])]["z_min"]), valves)
            map(x -> x["z_max"] = maximum(ts[string(x["index"])]["z_max"]), valves)
        end
    end
end


function set_valve_warm_start!(data::Dict{String, <:Any})
    apply_wm!(_set_valve_warm_start!, data)
end


function _set_valve_warm_start!(data::Dict{String, <:Any})
    for valve in values(data["valve"])
        flow_mid = 0.5 * (valve["flow_min"] + valve["flow_max"])

        valve["q_start"] = get(valve, "q", flow_mid)
        valve["qp_start"] = max(0.0, get(valve, "q", flow_mid))
        valve["qn_start"] = max(0.0, -get(valve, "q", flow_mid))

        valve["y_valve_start"] = get(valve, "q", 0.0) > 0.0 ? 1.0 : 0.0
        valve["z_valve_start"] = get(valve, "q", 0.0) > 0.0 ? 1.0 : 0.0
    end
end
