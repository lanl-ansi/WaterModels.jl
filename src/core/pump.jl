function _calc_pump_flow_min(pump::Dict{String,<:Any}, node_fr::Dict{String,Any}, node_to::Dict{String,Any})
    return max(0.0, get(pump, "flow_min", 0.0))
end


function _calc_pump_flow_min_forward(pump::Dict{String,<:Any}, node_fr::Dict{String,Any}, node_to::Dict{String,Any})
    flow_min_forward = get(pump, "flow_min_forward", _FLOW_MIN)
    return max(_calc_pump_flow_min(pump, node_fr, node_to), flow_min_forward)
end


function _calc_pump_flow_max_reverse(pump::Dict{String,<:Any}, node_fr::Dict{String,Any}, node_to::Dict{String,Any})
    flow_max_reverse = get(pump, "flow_max_reverse", 0.0)
    return min(_calc_pump_flow_max(pump, node_fr, node_to), flow_max_reverse)
end


function _calc_pump_flow_max(pump::Dict{String,<:Any}, node_fr::Dict{String,Any}, node_to::Dict{String,Any})
    coeff = _get_function_from_head_curve(pump["head_curve"])
    q_max_1 = (-coeff[2] + sqrt(coeff[2]^2 - 4.0*coeff[1]*coeff[3])) * inv(2.0*coeff[1])
    q_max_2 = (-coeff[2] - sqrt(coeff[2]^2 - 4.0*coeff[1]*coeff[3])) * inv(2.0*coeff[1])
    return min(max(q_max_1, q_max_2), get(pump, "flow_max", Inf))
end


function _calc_pump_flow_bounds_active(pump::Dict{String,<:Any})
    q_min, q_max = _calc_pump_flow_bounds(pump)
    q_min = max(max(get(pump, "flow_min_forward", _FLOW_MIN), _FLOW_MIN), q_min)
    return q_min, q_max
end
