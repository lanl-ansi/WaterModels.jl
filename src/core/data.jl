# Functions for working with the WaterModels internal data format.

function calc_head_bounds(junctions, reservoirs)
    junction_ids = [key for key in keys(junctions)]
    reservoir_ids = [key for key in keys(reservoirs)]
    ids = [junction_ids; reservoir_ids]

    head_min = Dict([(i, -Inf) for i in ids])
    head_max = Dict([(i, Inf) for i in ids])

    max_elev = maximum([junction["elev"] for junction in values(junctions)])
    max_head = maximum([reservoir["head"] for reservoir in values(reservoirs)])

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
        # Get the diameter of the pipe (meters).
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
    end

    return flow_min, flow_max
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
        diameter_vars = wm.var[:nw][wm.cnw][:psi][pipe_id]
        diameters = [key[1] for key in keys(diameter_vars)]
        for diameter in diameters
            if getvalue(wm.var[:nw][wm.cnw][:psi][pipe_id][diameter]) > 1.0e-4
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
