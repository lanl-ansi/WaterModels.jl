function aggregate_valves(subnetworks::Array{Dict{String, Any}, 1})
    return _aggregate_pipes(subnetworks, "valve")
end


function correct_valves!(data::Dict{String, <:Any})
    capacity = _calc_capacity_max(data)

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
end


function _calc_valve_flow_min(valve::Dict{String, <:Any}, capacity::Float64)
    flow_min_dir = valve["flow_direction"] == FLOW_DIRECTION_POSITIVE ? 0.0 : -Inf
    return max(-capacity, flow_min_dir, get(valve, "flow_min", -Inf))
end


function _calc_valve_flow_max(valve::Dict{String, <:Any}, capacity::Float64)
    flow_max_dir = valve["flow_direction"] == FLOW_DIRECTION_NEGATIVE ? 0.0 : Inf
    return min(capacity, flow_max_dir, get(valve, "flow_max", Inf))
end


function set_valve_warm_start!(data::Dict{String, <:Any})
    InfrastructureModels.apply!(_set_valve_warm_start!, data, wm_it_name)
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
