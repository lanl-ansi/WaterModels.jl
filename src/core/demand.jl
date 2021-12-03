function set_demand_bounds_from_time_series!(data::Dict{String,<:Any})
    wm_data = get_wm_data(data) # Get WaterModels-only portion of the data.
    @assert !ismultinetwork(wm_data) && haskey(wm_data, "time_series")
    func! = x -> _set_demand_bounds_from_time_series!(x, wm_data["time_series"])
    apply_wm!.(func!, values(wm_data["demand"]); apply_to_subnetworks = false)
end

function _set_demand_bounds_from_time_series!(demand::Dict{String,<:Any}, time_series::Dict{String,<:Any})
    # Get the index of the demand.
    demand_index = string(demand["index"])

    if haskey(time_series, "demand") && haskey(time_series["demand"], demand_index)
        # Get the time series data corresponding to the demand.
        demand_time_series = time_series["demand"][demand_index]
        
        if haskey(demand_time_series, "flow_min")
            # Set minimum flow to the minimum across all time.
            demand["flow_min"] = minimum(demand_time_series["flow_min"])
        end

        if haskey(demand_time_series, "flow_max")
            # Set the maximum flow to the maximum across all time.
            demand["flow_max"] = maximum(demand_time_series["flow_max"])
        end

        if haskey(demand_time_series, "flow_nominal")
            # Set the nominal flow to the mean across all time.
            demand["flow_nominal"] = Statistics.mean(demand_time_series["flow_nominal"])
        end

        # Ensure demand values are bounded as expected.
        @assert demand["flow_min"] <= demand["flow_nominal"]
        @assert demand["flow_nominal"] <= demand["flow_max"]
    end

    # If "flow_min" and "flow_max" are different, make the demand dispatchable.
    demand["dispatchable"] = demand["flow_min"] < demand["flow_max"]
end


function set_demand_warm_start!(data::Dict{String, <:Any})
    apply_wm!(_set_demand_warm_start!, data)
end


function _set_demand_warm_start!(data::Dict{String, <:Any})
    for demand in values(data["demand"])
        flow_mid = 0.5 * (demand["flow_min"] + demand["flow_max"])
        demand["q_demand_start"] = get(demand, "q", flow_mid)
    end
end
