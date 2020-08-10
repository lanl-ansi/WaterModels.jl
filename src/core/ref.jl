function _calc_head_loss_values(points::Array{Float64}, alpha::Float64)
    return [sign(x) * abs(x)^alpha for x in points]
end

function _calc_pump_gain_values(points::Array{Float64}, curve_fun::Array{Float64})
    return [curve_fun[1]*x*x + curve_fun[2]*x + curve_fun[3] for x in points]
end

function _calc_cubic_flow_values(points::Array{Float64}, curve_fun::Array{Float64})
    return [curve_fun[1]*x^3 + curve_fun[2]*x^2 + curve_fun[3]*x for x in points]
end

function _calc_efficiencies(points::Array{Float64}, curve::Array{Tuple{Float64,Float64}})
    q, eff = [[x[1] for x in curve], [x[2] for x in curve]]
    return Interpolations.LinearInterpolation(q, eff,
        extrapolation_bc=Interpolations.Flat()).(points)
end

function _get_function_from_head_curve(head_curve::Array{Tuple{Float64,Float64}})
    LsqFit.@. func(x, p) = p[1]*x*x + p[2]*x + p[3]

    if length(head_curve) > 1
        fit = LsqFit.curve_fit(func, first.(head_curve), last.(head_curve), [0.0, 0.0, 0.0])
        return LsqFit.coef(fit)
    elseif length(head_curve) == 1
        new_points = [(0.0, 1.33 * head_curve[1][2]), (2.0 * head_curve[1][1], 0.0)]
        head_curve = vcat(new_points, head_curve)
        fit = LsqFit.curve_fit(func, first.(head_curve), last.(head_curve), [0.0, 0.0, 0.0])
        return LsqFit.coef(fit)
    end
end


function calc_head_bounds(wm::AbstractWaterModel, n::Int=wm.cnw)
    # Compute the maximum elevation of all nodes in the network.
    max_head = maximum(node["elevation"] for (i, node) in ref(wm, n, :node))

    # Add contributions from tanks in the network to max_head.
    if length(ref(wm, n, :tank)) > 0
        max_head += sum(tank["max_level"] for (i, tank) in ref(wm, n, :tank))
    end

    # Add potential gains from pumps in the network to max_head.
    for (a, pump) in ref(wm, n, :pump)
        head_curve = ref(wm, n, :pump, a)["head_curve"]
        c = _get_function_from_head_curve(head_curve)
        max_head += c[3] - 0.25 * c[2]*c[2] * inv(c[1])
    end

    # Initialize the dictionaries for minimum and maximum heads.
    head_min = Dict((i, node["elevation"]) for (i, node) in ref(wm, n, :node))
    head_max = Dict((i, max_head) for (i, node) in ref(wm, n, :node))

    for (i, node) in ref(wm, n, :node)
        num_junctions = length(ref(wm, n, :node_junction, i))
        num_reservoirs = length(ref(wm, n, :node_reservoir, i))
        num_tanks = length(ref(wm, n, :node_tank, i))

        # If the node has zero demand, pressures can be negative.
        if num_junctions + num_reservoirs + num_tanks == 0
            head_min[i] -= 100.0
        end
    end

    for (i, reservoir) in ref(wm, n, :reservoir)
        # Fix head values at nodes with reservoirs to predefined values.
        node_id = reservoir["node"]
        head_min[node_id] = ref(wm, n, :node, node_id)["head"]
        head_max[node_id] = ref(wm, n, :node, node_id)["head"]
    end

    for (i, tank) in ref(wm, n, :tank)
        node_id = tank["node"]
        elevation = ref(wm, n, :node, node_id)["elevation"]
        head_min[node_id] = elevation + tank["min_level"]
        head_max[node_id] = elevation + tank["max_level"]
    end

    for (a, pressure_reducing_valve) in ref(wm, n, :pressure_reducing_valve)
        p_setting = ref(wm, n, :pressure_reducing_valve, a)["setting"]
        node_to = pressure_reducing_valve["node_to"]
        h_setting = ref(wm, n, :node, node_to)["elevation"] + p_setting
        head_min[node_to] = min(head_min[node_to], h_setting)
    end

    for (i, node) in ref(wm, n, :node)
        haskey(node, "h_min") && (head_min[i] = max(head_min[i], node["h_min"]))
        haskey(node, "h_max") && (head_max[i] = min(head_max[i], node["h_max"]))
    end

    # Return the dictionaries of lower and upper bounds.
    return head_min, head_max
end


function calc_flow_bounds(wm::AbstractWaterModel, n::Int=wm.cnw)
    h_lb, h_ub = calc_head_bounds(wm, n)
    alpha = ref(wm, n, :alpha)

    junctions = values(ref(wm, n, :junction))
    sum_demand = sum(abs(junction["demand"]) for junction in junctions)

    if :time_step in keys(ref(wm, n))
        time_step = ref(wm, n, :time_step)
        V_lb, V_ub = calc_tank_volume_bounds(wm, n)
        sum_demand += sum([(V_ub[i] - V_lb[i]) / time_step for i in ids(wm, n, :tank)])
    end

    lb, ub = Dict{String,Any}(), Dict{String,Any}()

    for name in ["check_valve", "pipe", "shutoff_valve"]
        lb[name], ub[name] = Dict{Int,Array{Float64}}(), Dict{Int,Array{Float64}}()

        for (a, comp) in ref(wm, n, Symbol(name))
            L, i, j = comp["length"], comp["node_fr"], comp["node_to"]
            resistances = ref(wm, n, :resistance, a)
            num_resistances = length(resistances)

            dh_lb, dh_ub = h_lb[i] - h_ub[j], h_ub[i] - h_lb[j]
            lb[name][a], ub[name][a] = zeros(num_resistances), zeros(num_resistances)

            for (r_id, r) in enumerate(resistances)
                lb[name][a][r_id] = sign(dh_lb) * (abs(dh_lb) * inv(L*r))^inv(alpha)
                lb[name][a][r_id] = max(lb[name][a][r_id], -sum_demand)

                ub[name][a][r_id] = sign(dh_ub) * (abs(dh_ub) * inv(L*r))^inv(alpha)
                ub[name][a][r_id] = min(ub[name][a][r_id], sum_demand)

                if comp["flow_direction"] == POSITIVE
                    lb[name][a][r_id] = max(lb[name][a][r_id], 0.0)
                elseif comp["flow_direction"] == NEGATIVE
                    ub[name][a][r_id] = min(ub[name][a][r_id], 0.0)
                end

                if haskey(comp, "diameters") && haskey(comp, "maximumVelocity")
                    D_a = comp["diameters"][r_id]["diameter"]
                    rate_bound = 0.25 * pi * comp["maximumVelocity"] * D_a * D_a
                    lb[name][a][r_id] = max(lb[name][a][r_id], -rate_bound)
                    ub[name][a][r_id] = min(ub[name][a][r_id], rate_bound)
                end
            end

            haskey(comp, "q_min") && (lb[name][a][1] = max(lb[name][a][1], comp["q_min"]))
            haskey(comp, "q_max") && (ub[name][a][1] = min(ub[name][a][1], comp["q_max"]))
        end
    end

    lb["pressure_reducing_valve"] = Dict{Int,Array{Float64}}()
    ub["pressure_reducing_valve"] = Dict{Int,Array{Float64}}()

    for (a, prv) in ref(wm, n, :pressure_reducing_valve)
        name = "pressure_reducing_valve"
        lb[name][a], ub[name][a] = [0.0], [sum_demand]
        haskey(prv, "q_min") && (lb[name][a][1] = max(lb[name][a][1], prv["q_min"]))
        haskey(prv, "q_max") && (ub[name][a][1] = min(ub[name][a][1], prv["q_max"]))
    end

    lb["pump"], ub["pump"] = Dict{Int,Array{Float64}}(), Dict{Int,Array{Float64}}()

    for (a, pump) in ref(wm, n, :pump)
        head_curve = ref(wm, n, :pump, a)["head_curve"]
        c = _get_function_from_head_curve(head_curve)
        q_max = (-c[2] + sqrt(c[2]^2 - 4.0*c[1]*c[3])) * inv(2.0*c[1])
        q_max = max(q_max, (-c[2] - sqrt(c[2]^2 - 4.0*c[1]*c[3])) * inv(2.0*c[1]))
        lb["pump"][a], ub["pump"][a] = [[0.0], [min(sum_demand, q_max)]]
        haskey(pump, "q_min") && (lb["pump"][a][1] = max(lb["pump"][a][1], pump["q_min"]))
        haskey(pump, "q_max") && (ub["pump"][a][1] = min(ub["pump"][a][1], pump["q_max"]))
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
