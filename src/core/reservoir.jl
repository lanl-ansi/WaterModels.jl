function set_reservoir_bounds_from_time_series!(data::Dict{String,<:Any})
    wm_data = get_wm_data(data) # Get WaterModels-only portion of the data.
    @assert !ismultinetwork(wm_data) && haskey(wm_data, "time_series")
    func! = x -> _set_reservoir_bounds_from_time_series!(x, wm_data["time_series"])
    apply_wm!.(func!, values(wm_data["reservoir"]); apply_to_subnetworks = false)
end


function _set_reservoir_bounds_from_time_series!(reservoir::Dict{String,<:Any}, time_series::Dict{String,<:Any})
    # Get the index of the reservoir.
    reservoir_index = string(reservoir["index"])

    if haskey(time_series, "reservoir") && haskey(time_series["reservoir"], reservoir_index)
        # Get the time series data corresponding to the reservoir.
        reservoir_time_series = time_series["reservoir"][reservoir_index]
 
        if haskey(reservoir_time_series, "head_min")
            # Set minimum head to the minimum across all time.
            reservoir["head_min"] = minimum(reservoir_time_series["head_min"])
        end

        if haskey(reservoir_time_series, "head_max")
            # Set the maximum head to the maximum across all time.
            reservoir["head_max"] = maximum(reservoir_time_series["head_max"])
        end

        if haskey(reservoir_time_series, "head_nominal")
            # Set the nominal head to the mean across all time.
            reservoir["head_nominal"] = Statistics.mean(reservoir_time_series["head_nominal"])
        end

        # Ensure reservoir values are bounded as expected.
        @assert reservoir["head_min"] <= reservoir["head_nominal"]
        @assert reservoir["head_nominal"] <= reservoir["head_max"]
    end

    # If "head_min" and "head_max" are different, make the reservoir dispatchable.
    reservoir["dispatchable"] = reservoir["head_min"] < reservoir["head_max"]
end


function set_reservoir_warm_start!(data::Dict{String, <:Any})
    apply_wm!(_set_reservoir_warm_start!, data)
end


function _set_reservoir_warm_start!(data::Dict{String, <:Any})
    for reservoir in values(data["reservoir"])
        reservoir["q_reservoir_start"] = get(reservoir, "q", 0.0)
    end
end
