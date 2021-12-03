function correct_regulators!(data::Dict{String,<:Any})
    apply_wm!(_correct_regulators!, data; apply_to_subnetworks = true)
end


function _correct_regulators!(data::Dict{String,<:Any})
    capacity = calc_capacity_max(data)

    for regulator in values(data["regulator"])
        _correct_status!(regulator)
        _correct_flow_direction!(regulator)
        _correct_regulator_flow_bounds!(regulator, capacity)
    end
end


function _correct_regulator_flow_bounds!(regulator::Dict{String,<:Any}, capacity::Float64)
    # Compute minimum flow for the regulator.
    flow_min_tmp = get(regulator, "flow_min", 0.0)
    regulator["flow_min"] = flow_min_tmp

    # Compute minimum _active_ (forward) flow for the regulator.
    flow_min_forward_tmp = get(regulator, "flow_min_forward", 0.0)
    regulator["flow_min_forward"] = max(flow_min_tmp, flow_min_forward_tmp)
    @assert regulator["flow_min"] <= regulator["flow_min_forward"]

    # Compute maximum flow for the regulator.
    flow_max_tmp = min(capacity, get(regulator, "flow_max", Inf))
    regulator["flow_max"] = flow_max_tmp
    @assert regulator["flow_min_forward"] <= regulator["flow_max"]
end