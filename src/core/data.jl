# Functions for working with the WaterModels internal data format.

function calc_flow_bounds(pipes)
    flow_min = Dict([(pipe_id, -Inf) for pipe_id in keys(pipes)])
    flow_max = Dict([(pipe_id, Inf) for pipe_id in keys(pipes)])

    for (pipe_id, pipe) in pipes
        # Diameter assumes original units of millimeters.
        diameter = pipe["diameter"] / 1000.0

        # A literature-based guess at the maximum velocity (meters per second).
        v_max = 10.0

        # Compute the flow bounds.
        max_absolute_flow = (pi / 4.0) * v_max * diameter^2
        flow_min[pipe_id] = -max_absolute_flow
        flow_max[pipe_id] = max_absolute_flow
    end

    return flow_min, flow_max
end

function calc_friction_factor(pipes, options)
    friction_factor = Dict([(pipe_id, -Inf) for pipe_id in keys(pipes)])
    headloss_type = options["headloss"]

    for (pipe_id, pipe) in pipes
        if headloss_type == "h-w"
            # Diameter assumes original units of millimeters.
            diameter = pipe["diameter"] / 1000.0

            # Roughness assumes no units.
            roughness = pipe["roughness"]

            # Length assumes original units of meters.
            length = pipe["length"]

            # Return the friction factor.
            friction_factor[pipe_id] = (10.67 * length) / (roughness^1.852 * diameter^4.87)
        elseif headloss_type == "d-w"
            # Diameter assumes original units of millimeters.
            diameter = pipe["diameter"] / 1000.0

            # Roughness assumes original units of millimeters.
            roughness = pipe["roughness"] / 1000.0

            # Length assumes original units of meters.
            length = pipe["length"]

            # Use standard gravitational acceleration on earth.
            g = 9.80665

            # Use the Prandtl-Kármán friction factor.
            f_s = 0.25 / log((roughness / diameter) / 3.71)^2

            # Return the friction factor.
            friction_factor[pipe_id] = (8.0 * length) / (pi^2 * g * diameter^5) * f_s
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
