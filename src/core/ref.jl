function calc_head_bounds(wm::GenericWaterModel, n::Int=wm.cnw)
    nodes = ref(wm, n, :nodes)
    tanks = ref(wm, n, :tanks)
    reservoirs = ref(wm, n, :reservoirs)

    if length(tanks) > 0 && length(reservoirs) > 0
        max_tank_elev = maximum(nodes[i]["elevation"] + tank["max_level"] for (i, tank) in tanks)
        max_reservoir_elev = maximum(res["elevation"] for (i, res) in reservoirs)
        max_elev = max_tank_elev + max_reservoir_elev
    elseif length(tanks) > 0 && length(reservoirs) == 0
        max_elev = maximum(nodes[i]["elevation"] + tank["max_level"] for (i, tank) in tanks)
    elseif length(tanks) == 0 && length(reservoirs) > 0
        max_elev = maximum(res["elevation"] for (i, res) in reservoirs)
    end

    # Initialize the dictionaries for minimum and maximum heads.
    head_min = Dict((i, -Inf) for (i, node) in nodes)
    head_max = Dict((i,  Inf) for (i, node) in nodes)

    for (i, node) in nodes
        # The minimum head at junctions must be above the initial elevation.
        if haskey(node, "minimumHead")
            head_min[i] = max(node["elevation"], node["minimumHead"])
        else
            head_min[i] = node["elevation"]
        end

        # The maximum head at junctions must be below the max reservoir height.
        if haskey(node, "maximumHead")
            head_max[i] = max(max_elev, node["maximumHead"])
        else
            # TODO: Is there a better general bound, here?
            head_max[i] = max_elev
        end
    end

    for (i, reservoir) in ref(wm, n, :reservoirs)
        # Head values at reservoirs are fixed.
        node_id = reservoir["reservoirs_node"]
        # TODO: Elevation should be a node attribute only.
        node = ref(wm, n, :reservoirs, node_id)
        head_min[node_id] = node["elevation"]
        head_max[node_id] = node["elevation"]
    end

    for (i, tank) in ref(wm, n, :tanks)
        node_id = tank["tanks_node"]
        node = ref(wm, n, :nodes, node_id)
        head_min[node_id] = node["elevation"] + tank["min_level"]
        head_max[node_id] = node["elevation"] + tank["max_level"]
    end

    # Return the dictionaries of lower and upper bounds.
    return head_min, head_max
end

function calc_head_difference_bounds(wm::GenericWaterModel, n::Int = wm.cnw)
    # Get placeholders for junctions and reservoirs.
    links = ref(wm, n, :links)

    # Initialize the dictionaries for minimum and maximum head differences.
    head_lbs, head_ubs = calc_head_bounds(wm, n)
    head_diff_min = Dict([(a, -Inf) for a in keys(links)])
    head_diff_max = Dict([(a, Inf) for a in keys(links)])

    # Compute the head difference bounds.
    for (a, link) in links
        head_diff_min[a] = head_lbs[link["f_id"]] - head_ubs[link["t_id"]]
        head_diff_max[a] = head_ubs[link["f_id"]] - head_lbs[link["t_id"]]
    end

    # Return the head difference bound dictionaries.
    return head_diff_min, head_diff_max
end

function calc_flow_rate_bounds(wm::GenericWaterModel, n::Int=wm.cnw)
    links = ref(wm, n, :links)
    dh_lb, dh_ub = calc_head_difference_bounds(wm, n)

    alpha = ref(wm, n, :alpha)
    junctions = values(ref(wm, n, :junctions))
    sum_demand = sum(junction["demand"] for junction in junctions)

    lb = Dict([(a, Float64[]) for a in keys(links)])
    ub = Dict([(a, Float64[]) for a in keys(links)])

    for (a, pipe) in ref(wm, n, :pipes)
        L = pipe["length"]
        resistances = ref(wm, n, :resistance, a)
        num_resistances = length(resistances)

        lb[a] = zeros(Float64, (num_resistances,))
        ub[a] = zeros(Float64, (num_resistances,))

        for (r_id, r) in enumerate(resistances)
            lb[a][r_id] = -100.0
            ub[a][r_id] = 100.0

            if uppercase(pipe["status"]) == "CV"
                lb[a][r_id] = 0.0
            end

            #lb[a][r_id] = sign(dh_lb[a]) * (abs(dh_lb[a]) / (L * r))^(inv(alpha))
            #ub[a][r_id] = sign(dh_ub[a]) * (abs(dh_ub[a]) / (L * r))^(inv(alpha))

            ## TODO: These seem to be valid bounds when tanks and pumps aren't
            ## present, but it should be fixed to include these, too.
            #lb[a][r_id] = max(lb[a][r_id], -sum_demand)
            #ub[a][r_id] = min(ub[a][r_id], sum_demand)

            #if pipe["flow_direction"] == POSITIVE
            #    lb[a][r_id] = max(lb[a][r_id], 0.0)
            #elseif pipe["flow_direction"] == NEGATIVE
            #    ub[a][r_id] = min(ub[a][r_id], 0.0)
            #end

            #if haskey(pipe, "diameters") && haskey(pipe, "maximumVelocity")
            #    D_a = pipe["diameters"][r_id]["diameter"]
            #    v_a = pipe["maximumVelocity"]
            #    rate_bound = 0.25 * pi * v_a * D_a * D_a
            #    lb[a][r_id] = max(lb[a][r_id], -rate_bound)
            #    ub[a][r_id] = min(ub[a][r_id], rate_bound)
            #end
        end
    end

    for (a, pump) in ref(wm, n, :pumps)
        # TODO: Need better bounds, here.
        lb[a] = [0.0]
        ub[a] = [Inf]
    end

    return lb, ub
end

function calc_directed_flow_upper_bounds(wm::GenericWaterModel, alpha::Float64, n::Int=wm.cnw)
    # Get a dictionary of resistance values.
    dh_lb, dh_ub = calc_head_difference_bounds(wm, n)

    links = ref(wm, n, :links)
    ub_n = Dict([(a, Float64[]) for a in keys(links)])
    ub_p = Dict([(a, Float64[]) for a in keys(links)])

    junctions = values(ref(wm, n, :junctions))
    sum_demand = sum(junction["demand"] for junction in junctions)

    for (a, link) in links
        L = link["length"]
        R_a = ref(wm, n, :resistance, a)

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

function calc_tank_volume_bounds(wm::GenericWaterModel, n::Int=wm.cnw)
    lb = Dict{Int64, Float64}(i => 0.0 for i in ids(wm, n, :tanks))
    ub = Dict{Int64, Float64}(i => Inf for i in ids(wm, n, :tanks))

    for (i, tank) in ref(wm, n, :tanks)
        if !("curve_name" in keys(tank))
            surface_area = 0.25 * pi * tank["diameter"]^2
            lb[i] = min(tank["min_vol"], surface_area * tank["min_level"])
            ub[i] = surface_area * tank["max_level"]
        else
            Memento.error(_LOGGER, "Only cylindrical tanks are currently supported.")
        end
    end

    return lb, ub
end
