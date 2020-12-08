function correct_regulators!(data::Dict{String, <:Any})
    capacity = _calc_capacity_max(data)

    for (idx, regulator) in data["regulator"]
        _correct_regulator_flow_bounds!(regulator, capacity)
    end
end


function _correct_regulator_flow_bounds!(regulator::Dict{String, <:Any}, capacity::Float64)
    flow_min = _calc_regulator_flow_min(regulator, capacity)
    flow_max = _calc_regulator_flow_max(regulator, capacity)
    regulator["flow_min"], regulator["flow_max"] = flow_min, flow_max
    regulator["flow_min_forward"] = max(flow_min, get(regulator, "flow_min_forward", 0.0))
    regulator["flow_max_reverse"] = min(flow_max, get(regulator, "flow_max_reverse", 0.0))
end


function _calc_regulator_flow_min(regulator::Dict{String, <:Any}, capacity::Float64)
    return max(-capacity, get(regulator, "flow_min", 0.0))
end


function _calc_regulator_flow_max(regulator::Dict{String, <:Any}, capacity::Float64)
    return min(capacity, get(regulator, "flow_max", Inf))
end
