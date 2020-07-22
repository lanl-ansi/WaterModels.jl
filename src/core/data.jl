# Functions for working with the WaterModels data format.
function calc_resistance_hw(diameter::Float64, roughness::Float64)
    return 7.8828 * inv(0.849^1.852 * roughness^1.852 * diameter^4.8704)
end


"""
Turns given single network data into multinetwork data with a `count` replicate
of the given network. Note that this function performs a deepcopy of the
network data. Significant multinetwork space savings can often be achieved by
building application specific methods of building multinetwork with minimal
data replication.
"""
function replicate(data::Dict{String,<:Any}, count::Int; global_keys::Set{String}=Set{String}())
    return _IM.replicate(data, count, union(global_keys, _wm_global_keys))
end


function calc_resistances_hw(pipes::Dict{<:Any, <:Any})
    resistances = Dict([(a, Array{Float64, 1}()) for a in keys(pipes)])

    for (a, pipe) in pipes
        if haskey(pipe, "resistances")
            resistances[a] = sort(pipe["resistances"], rev=true)
        elseif haskey(pipe, "resistance")
            resistances[a] = vcat(resistances[a], pipe["resistance"])
        elseif haskey(pipe, "diameters")
            for entry in pipe["diameters"]
                r = calc_resistance_hw(entry["diameter"], pipe["roughness"])
                resistances[a] = vcat(resistances[a], r)
            end

            resistances[a] = sort(resistances[a], rev=true)
        else
            r = calc_resistance_hw(pipe["diameter"], pipe["roughness"])
            resistances[a] = vcat(resistances[a], r)
        end
    end

    return resistances
end


function calc_resistance_dw(diameter::Float64, roughness::Float64, viscosity::Float64, speed::Float64, density::Float64)
    # Compute Reynold's number.
    reynolds_number = density * speed * diameter * inv(viscosity)

    # Use the same Colebrook formula as in EPANET.
    w = 0.25 * pi * reynolds_number
    y1 = 4.61841319859 * inv(w^0.9)
    y2 = (roughness * inv(diameter)) * inv(3.7 * diameter) + y1
    y3 = -8.685889638e-01 * log(y2)
    return 0.0826 * inv(diameter^5) * inv(y3*y3)
end


function calc_resistances_dw(pipes::Dict{<:Any, <:Any}, viscosity::Float64)
    resistances = Dict([(a, Array{Float64, 1}()) for a in keys(pipes)])

    for (a, pipe) in pipes
        if haskey(pipe, "resistances")
            resistances[a] = sort(pipe["resistances"], rev=true)
        elseif haskey(pipe, "resistance")
            resistance = pipe["resistance"]
            resistances[a] = vcat(resistances[a], resistance)
        elseif haskey(pipe, "diameters")
            for entry in pipe["diameters"]
                # Get relevant values to compute the friction factor.
                diameter = entry["diameter"]
                roughness = pipe["roughness"]
                r = calc_resistance_dw(diameter, roughness, viscosity, 10.0, 1000.0)
                resistances[a] = vcat(resistances[a], r)
            end

            resistances[a] = sort(resistances[a], rev = true)
        elseif haskey(pipe, "friction_factor")
            # Return the overall friction factor.
            diameter = pipe["diameter"]
            resistances[a] = [0.0826551*inv(diameter^5) * pipe["friction_factor"]]
        else
            # Get relevant values to compute the friction factor.
            diameter = pipe["diameter"]
            roughness = pipe["roughness"]
            r = calc_resistance_dw(diameter, roughness, viscosity, 10.0, 1000.0)
            resistances[a] = vcat(resistances[a], r)
        end
    end

    return resistances
end


function calc_resistance_costs_hw(pipes::Dict{Int, <:Any})
    # Create placeholder costs dictionary.
    costs = Dict([(a, Array{Float64, 1}()) for a in keys(pipes)])

    for (a, pipe) in pipes
        if haskey(pipe, "diameters")
            resistances = Array{Float64, 1}()

            for entry in pipe["diameters"]
                resistance = calc_resistance_hw(entry["diameter"], pipe["roughness"])
                resistances = vcat(resistances, resistance)
                costs[a] = vcat(costs[a], entry["costPerUnitLength"])
            end

            sort_indices = sortperm(resistances, rev=true)
            costs[a] = costs[a][sort_indices]
        else
            costs[a] = vcat(costs[a], 0.0)
        end
    end

    return costs
end


function calc_resistance_costs_dw(pipes::Dict{Int, <:Any}, viscosity::Float64)
    # Create placeholder costs dictionary.
    costs = Dict([(a, Array{Float64, 1}()) for a in keys(pipes)])

    for (a, pipe) in pipes
        if haskey(pipe, "diameters")
            resistances = Array{Float64, 1}()

            for entry in pipe["diameters"]
                diameter = entry["diameter"]
                roughness = pipe["roughness"]
                resistance = calc_resistance_dw(diameter, roughness, viscosity, 10.0, 1000.0)
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


function calc_resistances(pipes::Dict{<:Any, <:Any}, viscosity::Float64, head_loss_type::String)
    if head_loss_type == "H-W"
        return calc_resistances_hw(pipes)
    elseif head_loss_type == "D-W"
        return calc_resistances_dw(pipes, viscosity)
    else
        Memento.error(_LOGGER, "Head loss formulation type \"$(head_loss_type)\" is not recognized.")
    end
end


function calc_resistance_costs(pipes::Dict{Int, <:Any}, viscosity::Float64, head_loss_type::String)
    if head_loss_type == "H-W"
        return calc_resistance_costs_hw(pipes)
    elseif head_loss_type == "D-W"
        return calc_resistance_costs_dw(pipes, viscosity)
    else
        Memento.error(_LOGGER, "Head loss formulation type \"$(head_loss_type)\" is not recognized.")
    end
end


function has_known_flow_direction(comp::Pair{Int, <:Any})
    return comp.second["flow_direction"] != UNKNOWN
end


function is_des_pipe(pipe::Pair{Int64, <:Any})
    return any([x in ["diameters", "resistances"] for x in keys(pipe.second)])
end


function is_des_pipe(pipe::Pair{String, <:Any})
    return any([x in ["diameters", "resistances"] for x in keys(pipe.second)])
end


"turns a single network and a time_series data block into a multi-network"
function make_multinetwork(data::Dict{String, <:Any}; global_keys::Set{String}=Set{String}())
    return InfrastructureModels.make_multinetwork(data, union(global_keys, _wm_global_keys))
end


function set_start_head!(data)
    for (i, node) in data["node"]
        node["h_start"] = node["h"]
    end
end


function set_start_reservoir!(data)
    for (i, reservoir) in data["reservoir"]
        reservoir["qr_start"] = reservoir["qr"]
    end
end


function set_start_undirected_flow_rate!(data::Dict{String, <:Any})
    for (a, pipe) in data["pipe"]
        pipe["q_start"] = pipe["q"]
    end
end


function set_start_directed_flow_rate!(data::Dict{String, <:Any})
    for (a, pipe) in data["pipe"]
        pipe["qn_start"] = pipe["q"] < 0.0 ? abs(pipe["q"]) : 0.0
        pipe["qp_start"] = pipe["q"] >= 0.0 ? abs(pipe["q"]) : 0.0
    end
end


function set_start_directed_head_difference!(data::Dict{String, <:Any})
    for (a, pipe) in data["pipe"]
        i, j = data["pipe"][a]["node_fr"], data["pipe"][a]["node_to"]
        dh = data["node"][string(i)]["h"] - data["node"][string(j)]["h"]
        pipe["dhp_start"] = max(0.0, dh)
        pipe["dhn_start"] = max(0.0, -dh)
    end
end


function set_start_resistance_des!(data::Dict{String, <:Any})
    viscosity = data["viscosity"]
    head_loss_type = data["head_loss"]
    resistances = calc_resistances(data["pipe"], viscosity, head_loss_type)

    for (a, pipe) in filter(is_des_pipe, data["pipe"])
        num_resistances = length(resistances[a])
        pipe["x_res_start"] = zeros(Float64, num_resistances)
        r_id, val = findmax(pipe["x_res_start"])
        pipe["x_res_start"][r_id] = 1.0
    end
end


function set_start_undirected_flow_rate_des!(data::Dict{String, <:Any})
    viscosity = data["viscosity"]
    head_loss_type = data["head_loss"]
    resistances = calc_resistances(data["pipe"], viscosity, head_loss_type)

    for (a, pipe) in filter(is_des_pipe, data["pipe"])
        num_resistances = length(resistances[a])
        pipe["q_des_start"] = zeros(Float64, num_resistances)
        r_id, val = findmax(pipe["x_res_start"])
        pipe["q_des_start"][r_id] = pipe["q"]
    end
end


function set_start_directed_flow_rate_des!(data::Dict{String, <:Any})
    viscosity = data["viscosity"]
    head_loss_type = data["head_loss"]
    resistances = calc_resistances(data["pipe"], viscosity, head_loss_type)

    for (a, pipe) in filter(is_des_pipe, data["pipe"])
        num_resistances = length(resistances[a])
        pipe["qp_des_start"] = zeros(Float64, num_resistances)
        pipe["qn_des_start"] = zeros(Float64, num_resistances)

        r_id = findfirst(r -> isapprox(r, pipe["r"], rtol=1.0e-4), resistances[a])
        pipe["qp_des_start"][r_id] = pipe["q"] >= 0.0 ? abs(pipe["q"]) : 0.0
        pipe["qn_des_start"][r_id] = pipe["q"] < 0.0 ? abs(pipe["q"]) : 0.0
    end
end


function set_start_flow_direction!(data::Dict{String, <:Any})
    for (a, pipe) in data["pipe"]
        pipe["y_start"] = pipe["q"] >= 0.0 ? 1.0 : 0.0
    end
end


function set_start_all!(data::Dict{String, <:Any})
    set_start_head!(data)
    set_start_directed_head_difference!(data)
    set_start_reservoir!(data)
    set_start_resistance_des!(data)
    set_start_directed_flow_rate_des!(data)
    set_start_undirected_flow_rate_des!(data)
    set_start_flow_direction!(data)
end
