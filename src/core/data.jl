# Functions for working with the WaterModels internal data format.

function calc_head_bounds(wm::GenericWaterModel, n::Int = wm.cnw)
    # Get indices of nodes used in the network.
    junction_ids = collect(ids(wm, :junctions))
    reservoir_ids = collect(ids(wm, :reservoirs))
    nodes = [junction_ids; reservoir_ids]

    # Get placeholders for junctions and reservoirs.
    junctions = wm.ref[:nw][n][:junctions]
    reservoirs = wm.ref[:nw][n][:reservoirs]

    # Get maximum elevation/head values at nodes.
    max_elev = maximum([node["elev"] for node in values(junctions)])
    max_head = maximum([node["head"] for node in values(reservoirs)])

    # Initialize the dictionaries for minimum and maximum heads.
    head_min = Dict([(i, -Inf) for i in nodes])
    head_max = Dict([(i, Inf) for i in nodes])

    for (i, junction) in junctions
        # The minimum head at junctions must be above the initial elevation.
        if haskey(junction, "minimumHead")
            head_min[i] = max(junction["elev"], junction["minimumHead"])
        else
            head_min[i] = junction["elev"]
        end

        # The maximum head at junctions must be below the max reservoir height.
        if haskey(junction, "maximumHead")
            head_max[i] = max(max(max_elev, max_head), junction["maximumHead"])
        else
            head_max[i] = max(max_elev, max_head)
        end
    end

    for (i, reservoir) in reservoirs
        # Head values at reservoirs are fixed.
        head_min[i] = reservoir["head"]
        head_max[i] = reservoir["head"]
    end

    # Return the dictionaries of lower and upper bounds.
    return head_min, head_max
end

function calc_head_difference_bounds(wm::GenericWaterModel, n::Int = wm.cnw)
    # Get placeholders for junctions and reservoirs.
    links = wm.ref[:nw][n][:links]

    # Initialize the dictionaries for minimum and maximum head differences.
    head_lbs, head_ubs = calc_head_bounds(wm, n)
    head_diff_min = Dict([(a, -Inf) for a in keys(links)])
    head_diff_max = Dict([(a, Inf) for a in keys(links)])

    # Compute the head difference bounds.
    for (a, link) in links
        head_diff_min[a] = head_lbs[link["node1"]] - head_ubs[link["node2"]]
        head_diff_max[a] = head_ubs[link["node1"]] - head_lbs[link["node2"]]
    end

    # Return the head difference bound dictionaries.
    return head_diff_min, head_diff_max
end

function calc_directed_flow_upper_bounds(wm::GenericWaterModel, alpha::Float64=1.852, n::Int=wm.cnw)
    # Get a dictionary of resistance values.
    dh_lb, dh_ub = calc_head_difference_bounds(wm, n)

    links = wm.ref[:nw][n][:links]
    ub_n = Dict([(a, Float64[]) for a in keys(links)])
    ub_p = Dict([(a, Float64[]) for a in keys(links)])

    junctions = values(wm.ref[:nw][n][:junctions])
    sum_demand = sum(junction["demand"] for junction in junctions)

    for (a, link) in links
        L = link["length"]
        R_a = wm.ref[:nw][n][:resistance][a]

        ub_n[a] = zeros(Float64, (length(R_a),))
        ub_p[a] = zeros(Float64, (length(R_a),))

        for r in 1:length(R_a)
            ub_n[a][r] = abs(dh_lb[a] / (L * R_a[r]))^(1.0 / alpha)
            ub_n[a][r] = min(ub_n[a][r], sum_demand)

            ub_p[a][r] = abs(dh_ub[a] / (L * R_a[r]))^(1.0 / alpha)
            ub_p[a][r] = min(ub_p[a][r], sum_demand)

            if link["flow_direction"] == POSITIVE || dh_lb[a] >= 0.0
                ub_n[a][r] = 0.0
            elseif link["flow_direction"] == NEGATIVE || dh_ub[a] <= 0.0
                ub_p[a][r] = 0.0
            end

            if haskey(link, "diameters") && haskey(link, "maximumVelocity")
                D_a = link["diameters"][r]["diameter"]
                v_a = link["maximumVelocity"]
                rate_bound = 0.25 * pi * v_a * D_a * D_a
                ub_n[a][r] = min(ub_n[a][r], rate_bound)
                ub_p[a][r] = min(ub_p[a][r], rate_bound)
            end
        end
    end

    return ub_n, ub_p
end

function calc_resistance_per_length_hw(pipe)
    diameter = pipe["diameter"]
    roughness = pipe["roughness"]
    return 10.67 / (roughness^1.852 * diameter^4.87)
end

function get_node_ids(link::Dict{String, Any})
    i = link["node1"]
    j = link["node2"]
    return i, j
end

function calc_resistances_hw(wm::GenericWaterModel, n::Int = wm.cnw)
    # Get placeholders for junctions and reservoirs.
    links = wm.ref[:nw][n][:links]
    resistances = Dict([(a, Array{Float64, 1}()) for a in keys(links)])

    # Initialize the dictionaries for minimum and maximum head differences.
    for (a, link) in links
        if haskey(link, "diameters")
            for entry in link["diameters"]
                diameter = entry["diameter"]
                roughness = link["roughness"]
                r = 10.67 / (roughness^1.852 * diameter^4.87)
                resistances[a] = vcat(resistances[a], r)
            end

            resistances = sort(resistances, rev = true)
        else
            diameter = link["diameter"]
            roughness = link["roughness"]
            r = 10.67 / (roughness^1.852 * diameter^4.87)
            resistances[a] = vcat(resistances[a], r)
        end
    end

    return resistances
end

function calc_resistances_hw(links::Dict{Int, Any})
    # Get placeholders for junctions and reservoirs.
    resistances = Dict([(a, Array{Float64, 1}()) for a in keys(links)])

    # Initialize the dictionaries for minimum and maximum head differences.
    for (a, link) in links
        if haskey(link, "resistances")
            resistances[a] = sort(link["resistances"], rev = true)
        elseif haskey(link, "resistance")
            resistance = link["resistance"]
            resistances[a] = vcat(resistances[a], resistance)
        elseif haskey(link, "diameters")
            for entry in link["diameters"]
                diameter = entry["diameter"]
                roughness = link["roughness"]
                r = 10.67 / (roughness^1.852 * diameter^4.87)
                resistances[a] = vcat(resistances[a], r)
            end

            resistances[a] = sort(resistances[a], rev = true)
        else
            diameter = link["diameter"]
            roughness = link["roughness"]
            r = 10.67 / (roughness^1.852 * diameter^4.87)
            resistances[a] = vcat(resistances[a], r)
        end
    end

    return resistances
end

function calc_resistance_costs_hw(links::Dict{Int, Any})
    # Create placeholder costs dictionary.
    costs = Dict([(a, Array{Float64, 1}()) for a in keys(links)])

    # Initialize the dictionaries for minimum and maximum head differences.
    for (a, link) in links
        if haskey(link, "diameters")
            resistances = Array{Float64, 1}()

            for entry in link["diameters"]
                diameter = entry["diameter"]
                roughness = link["roughness"]
                resistance = 10.67 / (roughness^1.852 * diameter^4.87)
                resistances = vcat(resistances, resistance)
                costs[a] = vcat(costs[a], entry["costPerUnitLength"])
            end

            sort_indices = sortperm(resistances, rev = true)
            costs[a] = costs[a][sort_indices]
        else
            costs[a] = vcat(costs[a], 0.0)
        end
    end

    return costs
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

function get_maximum_diameter(link::Pair{Int, Any})
    if haskey(link.second, "diameters")
        entries = link.second["diameters"]
        return max([entry["diameter"] for entry in entries])
    else
        return link.second["diameter"]
    end
end


function has_known_flow_direction(link::Pair{Int, Any})
    return link.second["flow_direction"] != UNKNOWN
end

function is_ne_link(link::Pair{Int, Any})
    return haskey(link.second, "diameters") ||
           haskey(link.second, "resistances")
end

function is_out_node(i::Int)
    return function (link::Pair{Int, Any})
        return link.second["node1"] == i
    end
end

function is_in_node(i::Int)
    return function (link::Pair{Int, Any})
        return link.second["node2"] == i
    end
end
