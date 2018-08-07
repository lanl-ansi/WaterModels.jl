# Functions for working with the WaterModels internal data format.

function calc_flow_bounds(pipes)
    flow_min = Dict([(pipe_id, -Inf) for pipe_id in keys(pipes)])
    flow_max = Dict([(pipe_id, Inf) for pipe_id in keys(pipes)])

    for (pipe_id, pipe) in pipes
        # Get the diameter of the pipe (meters).
        diameter = pipe["diameter"]

        # A literature-based guess at the maximum velocity (meters per second).
        v_max = 10.0

        # Compute the flow bounds (cubic meters per second).
        max_absolute_flow = (pi / 4.0) * v_max * diameter^2
        flow_min[pipe_id] = -max_absolute_flow
        flow_max[pipe_id] = max_absolute_flow
    end

    return flow_min, flow_max
end

function calc_friction_factor(pipes, reynolds, options)
    friction_factor = Dict([(pipe_id, -Inf) for pipe_id in keys(pipes)])
    headloss_type = options["headloss"]

    for (pipe_id, pipe) in pipes
        diameter = pipe["diameter"]
        roughness = pipe["roughness"]
        length = pipe["length"]
        reynold = reynolds[pipe_id]

        if headloss_type == "h-w" # If Hazen-Williams...
            friction_factor[pipe_id] = (10.67 * length) / (roughness^1.852 * diameter^4.87)
        elseif headloss_type == "d-w" # If Darcy-Weisbach...
            # Use the same Colebrook formula as in EPANET.
            w = 0.25 * pi * reynold
            y1 = 4.61841319859 / w^0.9
            y2 = (roughness / diameter) / (3.7 * diameter) + y1
            y3 = -8.685889638e-01 * log(y2)
            f_s = 1.0 / y3^2

            # Compute the overall friction factor.
            # f_s = 0.25 / log((roughness / diameter) / 3.71 + 5.74 / reynold^0.9)^2
            friction_factor[pipe_id] = 0.0826 * length / diameter^5 * f_s
        else
            error("Could not find a valid \"headloss\" option type.")
        end
    end

    return friction_factor
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

function calc_reynolds_number(pipes, options)
    reynolds_number = Dict([(pipe_id, 0.0) for pipe_id in keys(pipes)])
    density = 1000.0 # Water density (kilograms per cubic meter).
    velocity = 10.0 # Estimate of velocity in the pipe (meters per second).
    viscosity = options["viscosity"] * 1.0e-3 # Viscosity is 10^(-3) Pascals * seconds.

    for (pipe_id, pipe) in pipes
        diameter = pipe["diameter"]
        reynolds_number[pipe_id] = density * velocity * diameter / viscosity
    end

    return reynolds_number
end

function update_flow_directions(data, wm) #solution)
    for (pipe_id, pipe) in data["pipes"]
        q = getvalue(wm.var[:nw][wm.cnw][:q][pipe_id])  # wm.solution["pipes"][pipe_id]["q"]
        pipe["flow_direction"] = q > 0.0 ? POSITIVE : NEGATIVE
    end
end

function reset_flow_directions(data)
    for (pipe_id, pipe) in data["pipes"]
        pipe["flow_direction"] = UNKNOWN
    end
end
