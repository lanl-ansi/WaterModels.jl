function _get_function_from_pump_curve(pump_curve::Array{Tuple{Float64,Float64}})
    LsqFit.@. func(x, p) = p[1]*x*x + p[2]*x + p[3]

    if length(pump_curve) > 1
        fit = LsqFit.curve_fit(func, first.(pump_curve), last.(pump_curve), [0.0, 0.0, 0.0])
        return LsqFit.coef(fit)
    elseif length(pump_curve) == 1
        new_points = [(0.0, 1.33 * pump_curve[1][2]), (2.0 * pump_curve[1][1], 0.0)]
        pump_curve = vcat(new_points, pump_curve)
        fit = LsqFit.curve_fit(func, first.(pump_curve), last.(pump_curve), [0.0, 0.0, 0.0])
        return LsqFit.coef(fit)
    end
end

function get_link_id(wm::AbstractWaterModel, i::Int, j::Int, n::Int=wm.cnw)
    arcs = vcat(ref(wm, n, :node_arc_fr, i), ref(wm, n, :node_arc_to, i))
    arc_id = findfirst(x -> x[2] in [i, j] && x[3] in [i, j], arcs)
    return arcs[arc_id][1]
end

function calc_head_bounds(wm::AbstractWaterModel, n::Int=wm.cnw)
    nodes = ref(wm, n, :node)
    tanks = ref(wm, n, :tank)
    reservoirs = ref(wm, n, :reservoir)

    if length(tanks) > 0
        max_elev = maximum(node["elevation"] for (i, node) in nodes)
        max_elev += maximum(tank["max_level"] for (i, tank) in tanks)
    else
        max_elev = maximum(node["elevation"] for (i, node) in nodes)
    end

    # Add head gains from pumps in the network to max_elev.
    for (a, pump) in ref(wm, n, :pump)
        pump_curve = ref(wm, n, :pump, a)["pump_curve"]
        c = _get_function_from_pump_curve(pump_curve)
        max_elev += c[3] - 0.25 * c[2]*c[2] * inv(c[1])
    end

    # Initialize the dictionaries for minimum and maximum heads.
    head_min = Dict((i, -Inf) for (i, node) in nodes)
    head_max = Dict((i, Inf) for (i, node) in nodes)

    for (i, node) in nodes
        # The minimum head at junctions must be above the initial elevation.
        if haskey(node, "minimumHead")
            head_min[i] = max(node["elevation"], node["minimumHead"])
        else
            head_min[i] = node["elevation"]
        end

        # The maximum head at junctions must be below the max elevation.
        if haskey(node, "maximumHead")
            head_max[i] = min(max_elev, node["maximumHead"])
        else
            # TODO: Is there a better general bound, here?
            head_max[i] = max_elev
        end
    end

    for (i, reservoir) in ref(wm, n, :reservoir)
        # Head values at reservoirs are fixed.
        node_id = reservoir["reservoir_node"]

        # TODO: Elevation should be a node attribute only.
        node = ref(wm, n, :reservoir, node_id)
        head_min[node_id] = node["elevation"]
        head_max[node_id] = node["elevation"]
    end

    for (i, tank) in ref(wm, n, :tank)
        node_id = tank["tank_node"]
        node = ref(wm, n, :node, node_id)
        head_min[node_id] = node["elevation"] + tank["min_level"]
        head_max[node_id] = node["elevation"] + tank["max_level"]
    end

    # Return the dictionaries of lower and upper bounds.
    return head_min, head_max
end

function calc_flow_bounds(wm::AbstractWaterModel, n::Int=wm.cnw)
    links = ref(wm, n, :link)
    h_lb, h_ub = calc_head_bounds(wm, n)
    alpha = ref(wm, n, :alpha)

    junctions = values(ref(wm, n, :junction))
    sum_demand = sum(abs(junction["demand"]) for junction in junctions)

    time_step_int = ref(wm, n, :option, "time")["hydraulic_timestep"]
    time_step = convert(Float64, time_step_int)
    V_lb, V_ub = calc_tank_volume_bounds(wm, n)
    sum_demand += sum([(V_ub[i] - V_lb[i]) / time_step for i in ids(wm, n, :tank)])

    lb = Dict([(a, Float64[]) for a in keys(links)])
    ub = Dict([(a, Float64[]) for a in keys(links)])

    for (a, pipe) in ref(wm, n, :pipe)
        L = pipe["length"]
        i, j = [pipe["node_fr"], pipe["node_to"]]
        dh_lb, dh_ub = [h_lb[i] - h_ub[j], h_ub[i] - h_lb[j]]
        resistances = ref(wm, n, :resistance, a)
        num_resistances = length(resistances)

        lb[a] = zeros(num_resistances)
        ub[a] = zeros(num_resistances)

        for (r_id, r) in enumerate(resistances)
            lb[a][r_id] = sign(dh_lb) * (abs(dh_lb) * inv(L*r))^inv(alpha)
            #lb[a][r_id] = max(lb[a][r_id], -sum_demand)

            ub[a][r_id] = sign(dh_ub) * (abs(dh_ub) * inv(L*r))^inv(alpha)
            #ub[a][r_id] = min(ub[a][r_id], sum_demand)

            if pipe["flow_direction"] == POSITIVE || has_check_valve(pipe)
                lb[a][r_id] = max(lb[a][r_id], 0.0)
            elseif pipe["flow_direction"] == NEGATIVE
                ub[a][r_id] = min(ub[a][r_id], 0.0)
            end

            if haskey(pipe, "diameters") && haskey(pipe, "maximumVelocity")
                D_a = pipe["diameters"][r_id]["diameter"]
                v_a = pipe["maximumVelocity"]
                rate_bound = 0.25 * pi * v_a * D_a * D_a
                lb[a][r_id] = max(lb[a][r_id], -rate_bound)
                ub[a][r_id] = min(ub[a][r_id], rate_bound)
            end
        end
    end

    for (a, pump) in ref(wm, n, :pump)
        pump_curve = ref(wm, n, :pump, a)["pump_curve"]
        c = _get_function_from_pump_curve(pump_curve)
        q_max = (-c[2] + sqrt(c[2]^2 - 4.0*c[1]*c[3])) * inv(2.0*c[1])
        q_max = max(q_max, (-c[2] - sqrt(c[2]^2 - 4.0*c[1]*c[3])) * inv(2.0*c[1]))
        lb[a], ub[a] = [[0.0], [q_max]]
    end

    return lb, ub
end

function calc_tank_volume_bounds(wm::AbstractWaterModel, n::Int=wm.cnw)
    lb = Dict{Int64, Float64}(i => 0.0 for i in ids(wm, n, :tank))
    ub = Dict{Int64, Float64}(i => Inf for i in ids(wm, n, :tank))

    for (i, tank) in ref(wm, n, :tank)
        if !("curve_name" in keys(tank))
            surface_area = 0.25 * pi * tank["diameter"]^2
            lb[i] = max(tank["min_vol"], surface_area * tank["min_level"])
            ub[i] = surface_area * tank["max_level"]
        else
            Memento.error(_LOGGER, "Only cylindrical tanks are currently supported.")
        end
    end

    return lb, ub
end
