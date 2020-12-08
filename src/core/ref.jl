"""
Builds the ref dictionary from the data dictionary. Additionally the ref dictionary would
contain fields populated by the optional vector of ref_extensions provided as a keyword
argument.
"""
function build_ref(
    data::Dict{String,<:Any};
    ref_extensions::Vector{<:Function} = Vector{Function}([]),
)
    return _IM.build_ref(
        data,
        ref_add_core!,
        _wm_global_keys,
        ref_extensions = ref_extensions,
    )
end


function _calc_pump_energy_points(wm::AbstractWaterModel, nw::Int, pump_id::Int, num_points::Int)
    pump = ref(wm, nw, :pump, pump_id)
    constant = _DENSITY * _GRAVITY * ref(wm, nw, :time_step)

    if get(wm.ext, :use_best_efficiency_form, false)
        curve_fun = _calc_pump_best_efficiency_head_gain_curve(pump)
    else
        curve_fun = _calc_pump_head_gain_curve(pump)
    end

    q_min, q_max = get(pump, "flow_min_forward", _FLOW_MIN), pump["flow_max"]
    q_build = range(q_min, stop = q_max, length = 100)
    f_build = _calc_cubic_flow_values(collect(q_build), curve_fun)

    if haskey(pump, "efficiency_curve")
        eff_curve = pump["efficiency_curve"]
        eff = _calc_efficiencies(collect(q_build), eff_curve)
    else
        eff = pump["efficiency"]
    end

    return q_build, constant .* inv.(eff) .* f_build
end


function _calc_pump_energy_ua(wm::AbstractWaterModel, nw::Int, pump_id::Int, q::Array{Float64, 1})
    q_true, f_true = _calc_pump_energy_points(wm, nw, pump_id, 100)
    f_interp = Interpolations.LinearInterpolation(q_true, f_true).(q)

    for i in 2:length(q)
        slope = (f_interp[i] - f_interp[i-1]) * inv(q[i] - q[i-1])
        true_ids = filter(x -> q_true[x] >= q[i-1] && q_true[x] <= q[i], 1:length(q_true))
        f_est_s = f_interp[i-1] .+ (slope .* (q_true[true_ids] .- q[i-1]))
        est_err = max(0.0, maximum(f_est_s .- f_true[true_ids]))
        f_interp[i-1:i] .-= est_err
    end

    return f_interp
end


function _calc_pump_energy(wm::AbstractWaterModel, nw::Int, pump_id::Int, q::Float64)
    pump = ref(wm, nw, :pump, pump_id)
    constant = _DENSITY * _GRAVITY * ref(wm, nw, :time_step)

    if get(wm.ext, :use_best_efficiency_form, false)
        pc = _calc_pump_best_efficiency_head_gain_curve(pump)
    else
        pc = _calc_pump_head_gain_curve(pump)
    end

    if haskey(pump, "efficiency_curve")
        eff_curve = pump["efficiency_curve"]
        eff = _calc_efficiencies([q], eff_curve)[1]
    else
        eff = pump["efficiency"]
    end

    return constant * inv(eff) * (pc[1]*q^3 + pc[2]*q^2 + pc[3]*q)
end


function _calc_head_loss_values(points::Array{Float64}, alpha::Float64)
    return [sign(x) * abs(x)^alpha for x in points]
end


function _calc_pump_gain_values(points::Array{Float64}, curve_fun::Array{Float64})
    return [curve_fun[1]*x*x + curve_fun[2]*x + curve_fun[3] for x in points]
end


function _calc_cubic_flow_values(points::Array{Float64}, curve_fun::Array{Float64})
    return [curve_fun[1]*x^3 + curve_fun[2]*x^2 + curve_fun[3]*x for x in points]
end


function _calc_efficiencies(points::Array{Float64}, curve::Array{Tuple{Float64, Float64}})
    q, eff = [[x[1] for x in curve], [x[2] for x in curve]]
    return Interpolations.LinearInterpolation(q, eff,
        extrapolation_bc=Interpolations.Flat()).(points)
end


function _get_function_from_head_curve(head_curve::Array{Tuple{Float64, Float64}})
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


function _get_maximum_head_gain(pump::Dict{String,<:Any}, use_best_efficiency::Bool)
    if use_best_efficiency
        c = _calc_pump_best_efficiency_head_gain_curve(pump)
    else
        c = _calc_pump_head_gain_curve(pump)
    end

    q_at_max = -c[2] * inv(2.0 * c[1]) > 0.0 ? -c[2] * inv(2.0 * c[1]) : 0.0
    return c[1] * q_at_max^2 + c[2] * q_at_max + c[3]
end


function calc_demand_bounds(wm::AbstractWaterModel, n::Int=wm.cnw)
    # Initialize the dictionaries for minimum and maximum heads.
    demand_min = Dict((i, -Inf) for (i, demand) in ref(wm, n, :dispatchable_demand))
    demand_max = Dict((i, Inf) for (i, demand) in ref(wm, n, :dispatchable_demand))

    for (i, demand) in ref(wm, n, :dispatchable_demand)
        demand_min[i] = max(demand_min[i], demand["demand_min"])
        demand_max[i] = min(demand_max[i], demand["demand_max"])
    end

    return demand_min, demand_max
end


function calc_head_bounds(wm::AbstractWaterModel, n::Int=wm.cnw)
    # Compute the maximum elevation of all nodes in the network.
    max_head = maximum(node["elevation"] for (i, node) in ref(wm, n, :node))

    # Add contributions from tanks in the network to max_head.
    if length(ref(wm, n, :tank)) > 0
        max_head += maximum(tank["max_level"] for (i, tank) in ref(wm, n, :tank))
    end

    if length(ref(wm, n, :pump)) > 0
        use_bep = get(wm.ext, :use_best_efficiency_form, false)
        max_head += maximum(_get_maximum_head_gain(pump, use_bep) for (i, pump) in ref(wm, n, :pump))
    end

    # Initialize the dictionaries for minimum and maximum heads.
    head_min = Dict((i, node["elevation"]) for (i, node) in ref(wm, n, :node))
    head_max = Dict((i, max_head) for (i, node) in ref(wm, n, :node))

    for (i, node) in ref(wm, n, :node)
        num_demands = length(ref(wm, n, :node_demand, i))
        num_reservoirs = length(ref(wm, n, :node_reservoir, i))
        num_tanks = length(ref(wm, n, :node_tank, i))

        # If the node has zero demand, pressures can be negative.
        if num_demands + num_reservoirs + num_tanks == 0
            head_min[i] -= 100.0
        end
    end

    for (i, reservoir) in ref(wm, n, :reservoir)
        # Fix head values at nodes with reservoirs to predefined values.
        node = ref(wm, n, :node, reservoir["node"])

        if haskey(node, "h_min") && haskey(node, "h_max")
            head_min[reservoir["node"]] = node["h_min"]
            head_max[reservoir["node"]] = node["h_max"]
        else
            head_min[reservoir["node"]] = node["head"]
            head_max[reservoir["node"]] = node["head"]
        end
    end

    for (i, tank) in ref(wm, n, :tank)
        node_id = tank["node"]
        elevation = ref(wm, n, :node, node_id)["elevation"]
        head_min[node_id] = elevation + tank["min_level"]
        head_max[node_id] = elevation + tank["max_level"]
    end

    for (a, regulator) in ref(wm, n, :regulator)
        p_setting = ref(wm, n, :regulator, a)["setting"]
        node_to = regulator["node_to"]
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

    nondispatchable_demands = values(ref(wm, n, :nondispatchable_demand))
    sum_demand = length(nondispatchable_demands) > 0 ? sum(demand["flow_nominal"] for demand in nondispatchable_demands) : 0.0
    dispatchable_demands = values(ref(wm, n, :dispatchable_demand))
    sum_demand += length(dispatchable_demands) > 0 ? sum(demand["flow_max"] for demand in dispatchable_demands) : 0.0

    if :time_step in keys(ref(wm, n))
        time_step = ref(wm, n, :time_step)
        V_lb, V_ub = calc_tank_volume_bounds(wm, n)
        sum_demand += sum([(V_ub[i] - V_lb[i]) / time_step for i in ids(wm, n, :tank)])
    end

    lb, ub = Dict{String,Any}(), Dict{String,Any}()

    for name in ["pipe", "des_pipe"]
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

            haskey(comp, "flow_min") && (lb[name][a][1] = max(lb[name][a][1], comp["flow_min"]))
            haskey(comp, "flow_max") && (ub[name][a][1] = min(ub[name][a][1], comp["flow_max"]))
        end
    end

    lb["regulator"] = Dict{Int,Array{Float64}}()
    ub["regulator"] = Dict{Int,Array{Float64}}()

    for (a, regulator) in ref(wm, n, :regulator)
        name = "regulator"
        lb[name][a], ub[name][a] = [0.0], [sum_demand]
        haskey(regulator, "flow_min") && (lb[name][a][1] = max(lb[name][a][1], regulator["flow_min"]))
        haskey(regulator, "flow_max") && (ub[name][a][1] = min(ub[name][a][1], regulator["flow_max"]))
    end

    lb["short_pipe"] = Dict{Int,Float64}()
    ub["short_pipe"] = Dict{Int,Float64}()

    for (a, short_pipe) in ref(wm, n, :short_pipe)
        lb["short_pipe"][a], ub["short_pipe"][a] = -sum_demand, sum_demand
        haskey(short_pipe, "flow_min") && (lb["short_pipe"][a] = max(lb["short_pipe"][a], short_pipe["flow_min"]))
        haskey(short_pipe, "flow_max") && (ub["short_pipe"][a] = min(ub["short_pipe"][a], short_pipe["flow_max"]))

        if short_pipe["flow_direction"] == POSITIVE
            lb["short_pipe"][a] = max(lb["short_pipe"][a], 0.0)
        elseif short_pipe["flow_direction"] == NEGATIVE
            ub["short_pipe"][a] = min(ub["short_pipe"][a], 0.0)
        end
    end

    lb["valve"] = Dict{Int,Float64}()
    ub["valve"] = Dict{Int,Float64}()

    for (a, valve) in ref(wm, n, :valve)
        lb["valve"][a], ub["valve"][a] = -sum_demand, sum_demand
        haskey(valve, "flow_min") && (lb["valve"][a] = max(lb["valve"][a], valve["flow_min"]))
        haskey(valve, "flow_max") && (ub["valve"][a] = min(ub["valve"][a], valve["flow_max"]))

        if valve["flow_direction"] == POSITIVE
            lb["valve"][a] = max(lb["valve"][a], 0.0)
        elseif valve["flow_direction"] == NEGATIVE
            ub["valve"][a] = min(ub["valve"][a], 0.0)
        end
    end

    lb["pump"], ub["pump"] = Dict{Int,Array{Float64}}(), Dict{Int,Array{Float64}}()

    for (a, pump) in ref(wm, n, :pump)
        if get(wm.ext, :use_best_efficiency_form, false)
            c = _calc_pump_best_efficiency_head_gain_curve(pump)
        else
            c = _calc_pump_head_gain_curve(pump)
        end

        q_max = (-c[2] + sqrt(c[2]^2 - 4.0*c[1]*c[3])) * inv(2.0*c[1])
        q_max = max(q_max, (-c[2] - sqrt(c[2]^2 - 4.0*c[1]*c[3])) * inv(2.0*c[1]))
        lb["pump"][a], ub["pump"][a] = [[0.0], [min(sum_demand, q_max)]]

        haskey(pump, "flow_min") && (lb["pump"][a][1] = max(lb["pump"][a][1], pump["flow_min"]))
        haskey(pump, "flow_max") && (ub["pump"][a][1] = min(ub["pump"][a][1], pump["flow_max"]))
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
