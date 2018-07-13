# Functions for working with the WaterModels internal data format.

function calc_flow_bounds(pipes, diameters)
    flow_min = Dict([(pipe_id, -Inf) for pipe_id in keys(pipes)])
    flow_max = Dict([(pipe_id, Inf) for pipe_id in keys(pipes)])

    for (pipe_id, pipe) in pipes
        # Get the diameter of the pipe (meters).
        diameter = diameters[pipe_id]

        # A literature-based guess at the maximum velocity (meters per second).
        v_max = 10.0

        # Compute the flow bounds (cubic meters per second).
        max_absolute_flow = (pi / 4.0) * v_max * diameter^2
        flow_min[pipe_id] = -max_absolute_flow
        flow_max[pipe_id] = max_absolute_flow
    end

    return flow_min, flow_max
end

function calc_demand(junctions, options)
    demand = Dict([(junction_id, 0.0) for junction_id in keys(junctions)])
    demand_units = options["units"]
    demand_multiplier = options["demand_multiplier"]

    for (junction_id, junction) in junctions
        if demand_units == "lps" # If liters per second...
            # Convert from liters per second to cubic meters per second.
            demand[junction_id] = junction["demand"] * 1.0e-3 * demand_multiplier
        elseif demand_units == "gpm" # If gallons per minute...
            # Convert from gallons per minute to cubic meters per second.
            demand[junction_id] = junction["demand"] * 6.30902e-5 * demand_multiplier
        else
            error("Could not find a valid \"units\" option type.")
        end
    end

    return demand
end

function calc_diameter(pipes, options)
    diameter = Dict([(pipe_id, 0.0) for pipe_id in keys(pipes)])
    demand_units = options["units"]

    for (pipe_id, pipe) in pipes
        if demand_units == "lps" # If liters per second...
            # Convert diameter from millimeters to meters.
            diameter[pipe_id] = 0.001 * pipe["diameter"]
        elseif demand_units == "gpm" # If gallons per minute...
            # Convert diameter from inches to meters.
            diameter[pipe_id] = 0.0254 * pipe["diameter"]
        else
            error("Could not find a valid \"units\" option type.")
        end
    end

    return diameter
end

function calc_elev(junctions, options)
    junction_ids = [key for key in keys(junctions)]
    elev = Dict([(junction_id, 0.0) for junction_id in junction_ids])
    demand_units = options["units"]

    for (junction_id, junction) in junctions
        if demand_units == "lps" # If liters per second...
            # Retain the original value (in meters).
            elev[junction_id] = junction["elev"]
        elseif demand_units == "gpm" # If gallons per minute...
            # Convert from feet to meters.
            elev[junction_id] = junction["elev"] * 0.3048
        else
            error("Could not find a valid \"units\" option type.")
        end
    end

    return elev
end

function calc_friction_factor(pipes, lengths, diameters, roughnesses, reynolds, options)
    friction_factor = Dict([(pipe_id, -Inf) for pipe_id in keys(pipes)])
    headloss_type = options["headloss"]

    for (pipe_id, pipe) in pipes
        diameter = diameters[pipe_id]
        roughness = roughnesses[pipe_id]
        length = lengths[pipe_id]
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

function calc_head(reservoirs, options)
    reservoir_ids = [key for key in keys(reservoirs)]
    head = Dict([(reservoir_id, 0.0) for reservoir_id in reservoir_ids])
    demand_units = options["units"]

    for (reservoir_id, reservoir) in reservoirs
        if demand_units == "lps" # If liters per second...
            # Retain the original value (in meters).
            head[reservoir_id] = reservoir["head"]
        elseif demand_units == "gpm" # If gallons per minute...
            # Convert from feet to meters.
            head[reservoir_id] = reservoir["head"] * 0.3048
        else
            error("Could not find a valid \"units\" option type.")
        end
    end

    return head
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

function calc_length(pipes, options)
    length = Dict([(pipe_id, 0.0) for pipe_id in keys(pipes)])
    demand_units = options["units"]

    for (pipe_id, pipe) in pipes
        if demand_units == "lps"
            # Retain the original value (in meters).
            length[pipe_id] = pipe["length"]
        elseif demand_units == "gpm" # If gallons per minute...
            # Convert length from feet to meters.
            length[pipe_id] = 0.3048 * pipe["length"]
        else
            error("Could not find a valid \"units\" option type.")
        end
    end

    return length
end

function calc_reynolds_number(pipes, diameters, options)
    reynolds_number = Dict([(pipe_id, 0.0) for pipe_id in keys(pipes)])
    density = 1000.0 # Water density (kilograms per cubic meter).
    velocity = 10.0 # Estimate of velocity in the pipe (meters per second).
    viscosity = options["viscosity"] * 1.0e-3 # Viscosity is 10^(-3) Pascals * seconds.

    for (pipe_id, pipe) in pipes
        diameter = diameters[pipe_id]
        reynolds_number[pipe_id] = density * velocity * diameter / viscosity
    end

    return reynolds_number
end

function calc_roughness(pipes, options)
    roughness = Dict([(pipe_id, 0.0) for pipe_id in keys(pipes)])
    headloss_type = options["headloss"]
    demand_units = options["units"]

    for (pipe_id, pipe) in pipes
        if headloss_type == "d-w" && demand_units == "gpm" 
            # Convert roughness from millifeet to meters.
            roughness[pipe_id] = 3.048e-4 * pipe["roughness"]
        elseif headloss_type == "d-w" && demand_units == "lps"
            # Convert roughness from millimeters to meters.
            roughness[pipe_id] = 0.001 * pipe["roughness"]
        elseif headloss_type == "h-w"
            # Retain the original value (unitless).
            roughness[pipe_id] = pipe["roughness"]
        else
            error("Could not find a valid \"headloss\" and \"units\" option combination.")
        end
    end

    return roughness
end
