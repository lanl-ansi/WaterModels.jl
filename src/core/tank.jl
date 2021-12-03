function set_tank_bounds_from_time_series!(data::Dict{String,<:Any})
    wm_data = get_wm_data(data) # Get WaterModels-only portion of the data.
    @assert !ismultinetwork(wm_data) && haskey(wm_data, "time_series")
    func! = x -> _set_tank_bounds_from_time_series!(x, wm_data["time_series"])
    apply_wm!.(func!, values(wm_data["tank"]); apply_to_subnetworks = false)
end


function _set_tank_bounds_from_time_series!(tank::Dict{String,<:Any}, time_series::Dict{String,<:Any})
    # Get the index of the tank.
    tank_index = string(tank["index"])

    if haskey(time_series, "tank") && haskey(time_series["tank"], tank_index)
        # Get the time series data corresponding to the tank.
        tank_time_series = time_series["tank"][tank_index]
        
        if haskey(tank_time_series, "min_level")
            # Set minimum level to the minimum across all time.
            tank["min_level"] = minimum(tank_time_series["min_level"])
        end

        if haskey(tank_time_series, "max_level")
            # Set the maximum level to the maximum across all time.
            tank["max_level"] = maximum(tank_time_series["max_level"])
        end

        # Ensure tank minimum and maximum levels are sensible.
        @assert tank["min_level"] <= tank["max_level"]
    end

    # If "min_level" and "max_level" are different, make the tank dispatchable.
    tank["dispatchable"] = tank["min_level"] < tank["max_level"]
end


function make_tank_start_dispatchable!(data::Dict{String,<:Any})
    wm_data = get_wm_data(data)

    if _IM.ismultinetwork(wm_data)
        nw_ids = sort(collect(keys(wm_data["nw"])))
        start_nw = string(sort([parse(Int, i) for i in nw_ids])[1])

        for tank in values(wm_data["nw"][start_nw]["tank"])
            tank["dispatchable"] = true
        end
    else
        for tank in values(wm_data["tank"])
            tank["dispatchable"] = true
        end
    end
end


function set_tank_warm_start!(data::Dict{String, <:Any})
    apply_wm!(_set_tank_warm_start!, data)
end


function _set_tank_warm_start!(data::Dict{String, <:Any})
    for tank in values(data["tank"])
        tank["q_tank_start"] = get(tank, "q", 0.0)
    end
end