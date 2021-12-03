function _set_link_bounds_from_time_series!(link::Dict{String,<:Any}, type::String, time_series::Dict{String,<:Any})
    # Get the index of the link.
    link_index = string(link["index"])

    if haskey(time_series, type) && haskey(time_series[type], link_index)
        # Get the time series data corresponding to the link.
        link_time_series = time_series[type][link_index]
        
        if haskey(link_time_series, "flow_min")
            # Set to the minimum flow across all time.
            link["flow_min"] = minimum(link_time_series["flow_min"])
        end

        if haskey(link_time_series, "flow_min_forward")
            # Set to the minimum positively-directed flow across all time.
            flow_min_forward = minimum(link_time_series["flow_min_forward"])
            link["flow_min_forward"] = flow_min_forward
            @assert link["flow_min"] <= link["flow_min_forward"]
        elseif haskey(link, "flow_min_forward")
            @assert link["flow_min"] <= link["flow_min_forward"]
        end

        if haskey(link_time_series, "flow_max")
            # Set to the maximum flow across all time.
            link["flow_max"] = maximum(link_time_series["flow_max"])
        end

        if haskey(link_time_series, "flow_max_reverse")
            # Set to the maximum negatively-directed flow across all time.
            link["flow_max_reverse"] = maximum(link_time_series["flow_max_reverse"])
            @assert link["flow_max_reverse"] <= link["flow_max"]
        end

        # Ensure link values are bounded as expected.
        @assert link["flow_min"] <= link["flow_max"]

        # Set flow direction based on bounds.
        if link["flow_max"] < 0.0
            link["flow_direction"] = FLOW_DIRECTION_NEGATIVE
        elseif link["flow_min"] > 0.0
            link["flow_direction"] = FLOW_DIRECTION_POSITIVE
        else
            link["flow_direction"] = FLOW_DIRECTION_UNKNOWN
        end
    end
end