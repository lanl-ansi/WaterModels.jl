# Functions for working with the WaterModels internal data format.

function calc_head_bounds(wm::GenericWaterModel, n::Int=wm.cnw)
    # Get indices of nodes used in the network.
    junction_ids = collect(ids(wm, :junctions))
    reservoir_ids = collect(ids(wm, :reservoirs))
    tank_ids = collect(ids(wm, :tanks))
    nodes = [junction_ids; reservoir_ids; tank_ids]

    # Get placeholders for junctions and reservoirs.
    junctions = ref(wm, n, :junctions)
    reservoirs = ref(wm, n, :reservoirs)
    tanks = ref(wm, n, :tanks)

    # Get maximum elevation/head values at nodes.
    max_elev = maximum([node["elevation"] for node in values(junctions)])
    max_head = maximum([node["head"] for node in values(reservoirs)])

    # Initialize the dictionaries for minimum and maximum heads.
    head_min = Dict([(i, -Inf) for i in nodes])
    head_max = Dict([(i, Inf) for i in nodes])

    for (i, junction) in junctions
        # The minimum head at junctions must be above the initial elevation.
        if haskey(junction, "minimumHead")
            head_min[i] = max(junction["elevation"], junction["minimumHead"])
        else
            head_min[i] = junction["elevation"]
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

    for (i, tank) in tanks
        head_min[i] = tank["elevation"] + tank["min_level"] 
        head_max[i] = tank["elevation"] + tank["max_level"] 
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
    pipes = ref(wm, n, :pipes)
    dh_lb, dh_ub = calc_head_difference_bounds(wm, n)

    alpha = ref(wm, n, :alpha)
    junctions = values(ref(wm, n, :junctions))
    sum_demand = sum(junction["demand"] for junction in junctions)

    lb = Dict([(a, Float64[]) for a in keys(links)])
    ub = Dict([(a, Float64[]) for a in keys(links)])

    for (a, link) in pipes
        L = link["length"]
        resistances = ref(wm, n, :resistance, a)
        num_resistances = length(resistances)

        lb[a] = zeros(Float64, (num_resistances,))
        ub[a] = zeros(Float64, (num_resistances,))

        for (r_id, r) in enumerate(resistances)
            lb[a][r_id] = sign(dh_lb[a]) * (abs(dh_lb[a]) / (L * r))^(inv(alpha))
            lb[a][r_id] = max(lb[a][r_id], -sum_demand)

            ub[a][r_id] = sign(dh_ub[a]) * (abs(dh_ub[a]) / (L * r))^(inv(alpha))
            ub[a][r_id] = min(ub[a][r_id], sum_demand)

            if link["flow_direction"] == POSITIVE
                lb[a][r_id] = max(lb[a][r_id], 0.0)
            elseif link["flow_direction"] == NEGATIVE
                ub[a][r_id] = min(ub[a][r_id], 0.0)
            end

            if haskey(link, "diameters") && haskey(link, "maximumVelocity")
                D_a = link["diameters"][r_id]["diameter"]
                v_a = link["maximumVelocity"]
                rate_bound = 0.25 * pi * v_a * D_a * D_a
                lb[a][r_id] = max(lb[a][r_id], -rate_bound)
                ub[a][r_id] = min(ub[a][r_id], rate_bound)
            end
        end
    end

    for (a, link) in ref(wm, n, :pumps)
        # TODO: Need better bounds here.
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

function calc_resistance_hw(diameter::Float64, roughness::Float64)
    return 10.67 * inv(roughness^1.852 * diameter^4.87)
end

function calc_resistances_hw(links::Dict{<:Any, Any})
    resistances = Dict([(a, Array{Float64, 1}()) for a in keys(links)])

    for (a, link) in links
        if haskey(link, "resistances")
            resistances[a] = sort(link["resistances"], rev=true)
        elseif haskey(link, "resistance")
            resistances[a] = vcat(resistances[a], link["resistance"])
        elseif haskey(link, "diameters")
            for entry in link["diameters"]
                r = calc_resistance_hw(entry["diameter"], link["roughness"])
                resistances[a] = vcat(resistances[a], r)
            end

            resistances[a] = sort(resistances[a], rev=true)
        else
            r = calc_resistance_hw(link["diameter"], link["roughness"])
            resistances[a] = vcat(resistances[a], r)
        end
    end

    return resistances
end

function get_num_resistances(link::Dict{String, Any})
    if haskey(link, "resistances")
        return length(link["resistances"])
    elseif haskey(link, "diameters")
        return length(link["diameters"])
    else
        return 1
    end
end

function calc_resistance_dw(length_::Float64, diameter::Float64, roughness::Float64, viscosity::Float64, speed::Float64, density::Float64)
    # Compute Reynold's number.
    reynolds_number = density * speed * diameter * inv(viscosity)

    # Use the same Colebrook formula as in EPANET.
    w = 0.25 * pi * reynolds_number
    y1 = 4.61841319859 * inv(w^0.9)
    y2 = (roughness * inv(diameter)) * inv(3.7 * diameter) + y1
    y3 = -8.685889638e-01 * log(y2)
    return 0.0826 * length_ * inv(diameter^5) * inv(y3*y3)
end

function calc_resistances_dw(links::Dict{<:Any, Any}, viscosity::Float64)
    resistances = Dict([(a, Array{Float64, 1}()) for a in keys(links)])

    for (a, link) in links
        length_ = link["length"]

        if haskey(link, "resistances")
            resistances[a] = sort(link["resistances"], rev = true)
        elseif haskey(link, "resistance")
            resistance = link["resistance"]
            resistances[a] = vcat(resistances[a], resistance)
        elseif haskey(link, "diameters")
            for entry in link["diameters"]
                # Get relevant values to compute the friction factor.
                diameter = entry["diameter"]
                roughness = link["roughness"]
                r = calc_resistance_dw(length_, diameter, roughness, viscosity, 10.0, 1000.0)
                resistances[a] = vcat(resistances[a], r)
            end

            resistances[a] = sort(resistances[a], rev = true)
        elseif haskey(link, "friction_factor")
            # Return the overall friction factor.
            resistances[a] = [0.0826 * length_ * inv(diameter^5) * pipe["friction_factor"]]
        else
            # Get relevant values to compute the friction factor.
            diameter = link["diameter"]
            roughness = link["roughness"]
            r = calc_resistance_dw(length_, diameter, roughness, viscosity, 10.0, 1000.0)
            resistances[a] = vcat(resistances[a], r)
        end
    end

    return resistances
end

function calc_resistance_costs_hw(links::Dict{Int, Any})
    # Create placeholder costs dictionary.
    costs = Dict([(a, Array{Float64, 1}()) for a in keys(links)])

    for (a, link) in links
        if haskey(link, "diameters")
            resistances = Array{Float64, 1}()

            for entry in link["diameters"]
                resistance = calc_resistance_hw(entry["diameter"], link["roughness"])
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

function calc_resistance_costs_dw(links::Dict{Int, Any}, viscosity::Float64)
    # Create placeholder costs dictionary.
    costs = Dict([(a, Array{Float64, 1}()) for a in keys(links)])

    for (a, link) in links
        length_ = link["length"]

        if haskey(link, "diameters")
            resistances = Array{Float64, 1}()

            for entry in link["diameters"]
                diameter = entry["diameter"]
                roughness = link["roughness"]
                resistance = calc_resistance_dw(length_, diameter, roughness, viscosity, 10.0, 1000.0)
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

function calc_resistances(links::Dict{<:Any, Any}, viscosity::Float64, head_loss_type::String)
    if head_loss_type == "H-W"
        return calc_resistances_hw(links)
    elseif head_loss_type == "D-W"
        return calc_resistances_dw(links, viscosity)
    else
        Memento.error(_LOGGER, "Head loss formulation type \"$(head_loss_type)\" is not recognized.")
    end
end

function calc_resistance_costs(links::Dict{Int, Any}, viscosity::Float64, head_loss_type::String)
    if head_loss_type == "H-W"
        return calc_resistance_costs_hw(links)
    elseif head_loss_type == "D-W"
        return calc_resistance_costs_dw(links, viscosity)
    else
        Memento.error(_LOGGER, "Head loss formulation type \"$(head_loss_type)\" is not recognized.")
    end
end

function has_known_flow_direction(link::Pair{Int, Any})
    return link.second["flow_direction"] != UNKNOWN
end

function is_ne_link(link::Pair{Int64, Any})
    return any([x in ["diameters", "resistances"] for x in keys(link.second)])
end

function is_ne_link(link::Pair{String, Any})
    return any([x in ["diameters", "resistances"] for x in keys(link.second)])
end

function is_out_node(i::Int)
    return function (link::Pair{Int, Any})
        return link.second["f_id"] == i
    end
end

function is_in_node(i::Int)
    return function (link::Pair{Int, Any})
        return link.second["t_id"] == i
    end
end


_wm_global_keys = Set(["time_series","mass_units","per_unit","options","flow_units"])

"""
Turns in given single network data in multinetwork data with a `count`
replicate of the given network. Note that this function performs a deepcopy of
the network data. Significant multinetwork space savings can often be achieved
by building application specific methods of building multinetwork with minimal
data replication.
"""
function replicate(sn_data::Dict{String,<:Any}, count::Int; global_keys::Set{String}=Set{String}())
    wm_global_keys = Set(["per_unit"])
    return InfrastructureModels.replicate(sn_data, count, global_keys=union(global_keys, _wm_global_keys))
end


"turns a single network and a time_series data block into a multi-network"
function make_multinetwork(data::Dict{String,<:Any})
    if InfrastructureModels.ismultinetwork(data)
        Memento.error(_LOGGER, "make_multinetwork does not support multinetwork data")
    end

    if !haskey(data, "time_series")
        Memento.error(_LOGGER, "make_multinetwork requires time_series data")
    end

    # hard coded for the moment
    steps = 24

    mn_data = InfrastructureModels.replicate(data, steps, global_keys=_wm_global_keys)
    time_series = pop!(mn_data, "time_series")

    for i in 1:steps
        nw_data = mn_data["nw"]["$(i)"]
        for (k,v) in data["time_series"]
            if isa(v, Dict) && haskey(nw_data, k)
                _update_data_timepoint!(nw_data[k], v, i)
            end
        end
    end

    return mn_data
end


"loads a single time point from a time_series data block into the current network"
function load_timepoint!(data::Dict{String,<:Any}, step::Int)
    if InfrastructureModels.ismultinetwork(data)
        Memento.error(_LOGGER, "load_timepoint! does not support multinetwork data")
    end

    if !haskey(data, "time_series")
        Memento.error(_LOGGER, "load_timepoint! requires time_series data")
    end

    for (k,v) in data["time_series"]
        if isa(v, Dict) && haskey(data, k)
            _update_data_timepoint!(data[k], v, step)
        end
    end

    data["step_index"] = step
end


"recursive call of _update_data"
function _update_data_timepoint!(data::Dict{String,<:Any}, new_data::Dict{String,<:Any}, step::Int)
    for (key, new_v) in new_data
        if haskey(data, key)
            v = data[key]
            if isa(v, Dict) && isa(new_v, Dict)
                _update_data_timepoint!(v, new_v, step)
            elseif (!isa(v, Dict) || !isa(v, Array)) && isa(new_v, Array)
                data[key] = new_v[step]
            else
                Memento.warn(_LOGGER, "skipping key $(key) because object types do not match, target $(typeof(v)) source $(typeof(new_v))")
            end
        else
            Memento.warn(_LOGGER, "skipping time_series key $(key) because it does not occur in the target data")
        end
    end
end




function set_start_head!(data)
    for (i, junction) in data["junctions"]
        junction["h_start"] = junction["h"]
    end
end

function set_start_undirected_flow_rate!(data::Dict{String, Any})
    for (a, pipe) in data["pipes"]
        pipe["q_start"] = pipe["q"]
    end
end

function set_start_directed_flow_rate!(data::Dict{String, Any})
    for (a, pipe) in data["pipes"]
        pipe["qn_start"] = pipe["q"] < 0.0 ? abs(pipe["q"]) : 0.0
        pipe["qp_start"] = pipe["q"] >= 0.0 ? abs(pipe["q"]) : 0.0
    end
end

function set_start_directed_head_difference!(data::Dict{String, Any})
    head_loss_type = data["options"]["headloss"]
    alpha = head_loss_type == "H-W" ? 1.852 : 2.0

    for (a, pipe) in data["pipes"]
        dh_abs = pipe["length"] * pipe["r"] * abs(pipe["q"])^(alpha)
        pipe["dhn_start"] = pipe["q"] < 0.0 ? dh_abs : 0.0
        pipe["dhp_start"] = pipe["q"] >= 0.0 ? dh_abs : 0.0
    end
end

function set_start_resistance_ne!(data::Dict{String, Any})
    viscosity = data["options"]["hydraulic"]["viscosity"]
    head_loss_type = data["options"]["hydraulic"]["headloss"]
    resistances = calc_resistances(data["pipes"], viscosity, head_loss_type)

    for (a, pipe) in filter(is_ne_link, data["pipes"])
        num_resistances = length(resistances[a])
        pipe["x_res_start"] = zeros(Float64, num_resistances)
        r_id = findfirst(r -> isapprox(r, pipe["r"], rtol=1.0e-4), resistances[a])
        pipe["x_res_start"][r_id] = 1.0
    end
end

function set_start_undirected_flow_rate_ne!(data::Dict{String, Any})
    viscosity = data["options"]["hydraulic"]["viscosity"]
    head_loss_type = data["options"]["hydraulic"]["headloss"]
    resistances = calc_resistances(data["pipes"], viscosity, head_loss_type)

    for (a, pipe) in filter(is_ne_link, data["pipes"])
        num_resistances = length(resistances[a])
        pipe["q_ne_start"] = zeros(Float64, num_resistances)
        r_id = findfirst(r -> isapprox(r, pipe["r"], rtol=1.0e-4), resistances[a])
        pipe["q_ne_start"][r_id] = pipe["q"]
    end
end

function set_start_directed_flow_rate_ne!(data::Dict{String, Any})
    viscosity = data["options"]["hydraulic"]["viscosity"]
    head_loss_type = data["options"]["hydraulic"]["headloss"]
    resistances = calc_resistances(data["pipes"], viscosity, head_loss_type)

    for (a, pipe) in filter(is_ne_link, data["pipes"])
        num_resistances = length(resistances[a])
        pipe["qn_ne_start"] = zeros(Float64, num_resistances)
        pipe["qp_ne_start"] = zeros(Float64, num_resistances)

        r_id = findfirst(r -> isapprox(r, pipe["r"], rtol=1.0e-4), resistances[a])
        pipe["qn_ne_start"][r_id] = pipe["q"] < 0.0 ? abs(pipe["q"]) : 0.0
        pipe["qp_ne_start"][r_id] = pipe["q"] >= 0.0 ? abs(pipe["q"]) : 0.0
    end
end

function set_start_flow_direction!(data::Dict{String, Any})
    for (a, pipe) in data["pipes"]
        pipe["x_dir_start"] = pipe["q"] >= 0.0 ? 1.0 : 0.0
    end
end
