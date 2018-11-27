# Functions for working with the WaterModels internal data format.

function calc_head_bounds(junctions, reservoirs)
    junction_ids = [key for key in keys(junctions)]
    reservoir_ids = [key for key in keys(reservoirs)]
    ids = [junction_ids; reservoir_ids]

    head_min = Dict([(i, -Inf) for i in ids])
    head_max = Dict([(i, Inf) for i in ids])

    max_elev = maximum([junc["elev"] for junc in values(junctions)])
    max_head = maximum([res["head"] for res in values(reservoirs)])

    for (i, junction) in junctions
        if haskey(junction, "minimumHead")
            head_min[i] = max(junction["elev"], junction["minimumHead"])
        else
            head_min[i] = junction["elev"]
        end

        if haskey(junction, "maximumHead")
            head_max[i] = max(max(max_elev, max_head), junction["maximumHead"])
        else
            head_max[i] = max(max_elev, max_head)
        end
    end

    for (i, reservoir) in reservoirs
        head_min[i] = reservoir["head"]
        head_max[i] = reservoir["head"]
    end

    return head_min, head_max
end

function calc_flow_bounds(pipes)
    flow_min = Dict([(pipe_id, -Inf) for pipe_id in keys(pipes)])
    flow_max = Dict([(pipe_id, Inf) for pipe_id in keys(pipes)])

    for (pipe_id, pipe) in pipes
        diameter = pipe["diameter"]
        max_absolute_flow = (pi / 4.0) * 10.0 * diameter^2

        if haskey(pipe, "minimumFlow")
            flow_min[pipe_id] = pipe["minimumFlow"]
        else
            flow_min[pipe_id] = -max_absolute_flow
        end

        if haskey(pipe, "maximumFlow")
            flow_max[pipe_id] = pipe["maximumFlow"]
        else
            flow_max[pipe_id] = max_absolute_flow
        end

        if pipe["flow_direction"] == POSITIVE
            flow_min[pipe_id] = max(0.0, flow_min[pipe_id])
        elseif pipe["flow_direction"] == NEGATIVE
            flow_max[pipe_id] = min(0.0, flow_max[pipe_id])
        end
    end

    return flow_min, flow_max
end

function calc_resistance_per_length_hw(pipe)
    diameter = pipe["diameter"]
    roughness = pipe["roughness"]
    return 10.67 / (roughness^1.852 * diameter^4.87)
end

function calc_friction_factor_hw(pipe)
    if haskey(pipe, "friction_factor")
        return pipe["friction_factor"]
    else
        diameter = pipe["diameter"]
        roughness = pipe["roughness"]
        length = pipe["length"]
        return (10.67 * length) / (roughness^1.852 * diameter^4.87)
    end
end

function calc_friction_factor_hw_ne(pipe, diameter)
    if haskey(pipe, "friction_factor")
        return pipe["friction_factor"]
    else
        roughness = pipe["roughness"]
        length = pipe["length"]
        return (10.67 * length) / (roughness^1.852 * diameter^4.87)
    end
end

function calc_friction_factor_dw(pipe, viscosity)
    diameter = pipe["diameter"]
    length = pipe["length"]

    if haskey(pipe, "friction_factor")
        # Return the overall friction factor.
        return 0.0826 * length / diameter^5 * pipe["friction_factor"]
    else
        # Get relevant values to compute the friction factor.
        roughness = pipe["roughness"]

        # Compute Reynold's number.
        density = 1000.0 # Water density (kilograms per cubic meter).
        velocity = 10.0 # Estimate of velocity in the pipe (meters per second).
        reynolds_number = density * velocity * diameter / viscosity

        # Use the same Colebrook formula as in EPANET.
        w = 0.25 * pi * reynolds_number
        y1 = 4.61841319859 / w^0.9
        y2 = (roughness / diameter) / (3.7 * diameter) + y1
        y3 = -8.685889638e-01 * log(y2)
        f_s = 1.0 / y3^2

        # Return the overall friction factor.
        return 0.0826 * length / diameter^5 * f_s
    end
end

function calc_head_difference_bounds(pipes)
    diff_min = Dict([(pipe_id, -Inf) for pipe_id in keys(pipes)])
    diff_max = Dict([(pipe_id, Inf) for pipe_id in keys(pipes)])

    for (pipe_id, pipe) in pipes
        # Compute the flow bounds.
        # TODO: Replace these with better bounds.
        diff_min[pipe_id] = -1000.0
        diff_max[pipe_id] = 1000.0
    end

    return diff_min, diff_max
end

function update_flow_directions(data, wm)
    for (pipe_id, pipe) in data["pipes"]
        q = getvalue(wm.var[:nw][wm.cnw][:q][pipe_id])
        pipe["flow_direction"] = q >= 0.0 ? POSITIVE : NEGATIVE
    end
end

function update_diameters(data, wm)
    for (pipe_id, pipe) in data["pipes"]
        a = parse(Int, pipe_id)
        entries = wm.data["pipes"][pipe_id]["diameters"]
        diameters = [entry["diameter"] for entry in entries]

        for diameter in diameters
            if getvalue(wm.var[:nw][wm.cnw][:psi][a][diameter]) > 1.0e-4
                pipe["diameter"] = diameter
                break
            end
        end
    end
end

function reset_flow_directions(data)
    for (pipe_id, pipe) in data["pipes"]
        pipe["flow_direction"] = UNKNOWN
    end
end

function set_maximum_diameter(pipe::Dict{String, Any})
    if haskey(pipe, "diameters")
        entries = pipe["diameters"]
        diameters = [entry["diameter"] for entry in entries]
        pipe["diameter"] = maximum(diameters)
        delete!(pipe, "diameters")
    end
end

function set_maximum_diameters(data::Dict{String, Any})
    for (a, pipe) in data["pipes"]
        set_maximum_diameter(pipe)
    end
end

function get_maximum_diameter(connection::Pair{Int, Any})
    if haskey(connection.second, "diameters")
        entries = connection.second["diameters"]
        return max([entry["diameter"] for entry in entries])
    else
        return connection.second["diameter"]
    end
end


function has_known_flow_direction(connection::Pair{Int, Any})
    return connection.second["flow_direction"] != UNKNOWN
end

function is_ne_pipe(connection::Pair{Int, Any})
    return haskey(connection.second, "diameters")
end

function is_out_node_function(i::Int)
    function is_out_node(connection::Pair{Int, Any})
        return parse(Int, connection.second["node1"]) == i
    end

    return is_out_node
end

function is_in_node_function(i::Int)
    function is_in_node(connection::Pair{Int, Any})
        return parse(Int, connection.second["node2"]) == i
    end

    return is_in_node
end
