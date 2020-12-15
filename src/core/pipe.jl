function correct_pipes!(data::Dict{String, <:Any})
    head_loss_form, visc = data["head_loss"], data["viscosity"]
    capacity = _calc_capacity_max(data)

    for (idx, pipe) in data["pipe"]
        # Get common connecting node data for later use.
        node_fr = data["node"][string(pipe["node_fr"])]
        node_to = data["node"][string(pipe["node_to"])]

        # Correct various pipe properties. The sequence is important, here.
        _correct_pipe_flow_bounds!(pipe, node_fr, node_to, head_loss_form, visc, capacity)
    end
end


function _correct_pipe_flow_bounds!(
    pipe::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any},
    form::String, viscosity::Float64, capacity::Float64)
    # Calculate flow bounds, which depend on the head loss formulation.
    flow_min = _calc_pipe_flow_min(pipe, node_fr, node_to, form, viscosity, capacity)
    flow_max = _calc_pipe_flow_max(pipe, node_fr, node_to, form, viscosity, capacity)
    pipe["flow_min"], pipe["flow_max"] = flow_min, flow_max
    pipe["flow_min_forward"] = max(flow_min, get(pipe, "flow_min_forward", 0.0))
    pipe["flow_max_reverse"] = min(flow_max, get(pipe, "flow_max_reverse", 0.0))
end


function _calc_pipe_flow_min(
    pipe::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any},
    form::String, viscosity::Float64, capacity::Float64)
    # Calculate minimum flow bound depending on the head loss formulation.
    if uppercase(form) == "H-W"
        return _calc_pipe_flow_min_hw(pipe, node_fr, node_to, capacity)
    elseif uppercase(form) == "D-W"
        return _calc_pipe_flow_min_dw(pipe, node_fr, node_to, viscosity, capacity)
    end
end


function _calc_pipe_flow_max(
    pipe::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any},
    form::String, viscosity::Float64, capacity::Float64)
    # Calculate maximum flow bound depending on the head loss formulation.
    if uppercase(form) == "H-W"
        return _calc_pipe_flow_max_hw(pipe, node_fr, node_to, capacity)
    elseif uppercase(form) == "D-W"
        return _calc_pipe_flow_max_dw(pipe, node_fr, node_to, viscosity, capacity)
    end
end


function _calc_pipe_flow_min_hw(
    pipe::Dict{String, <:Any}, node_fr::Dict{String, <:Any},
    node_to::Dict{String, <:Any}, capacity::Float64)
    # Calculate minimum flow bound per the Hazen-Williams head loss formulation.
    resistance = _calc_pipe_resistance_hw(pipe["diameter"], pipe["roughness"])
    loss = get(node_fr, "head_min", -Inf) - get(node_to, "head_max", Inf)
    flow_min_loss = sign(loss) * (abs(loss) * inv(pipe["length"] * resistance))^inv(1.852)
    flow_min_dir = pipe["flow_direction"] == POSITIVE ? 0.0 : -Inf
    return max(-capacity, flow_min_loss, flow_min_dir, get(pipe, "flow_min", -Inf))
end


function _calc_pipe_flow_min_dw(
    pipe::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any},
    viscosity::Float64, capacity::Float64)
    # Calculate minimum flow bound per the Hazen-Williams head loss formulation.
    diameter, roughness, speed = pipe["diameter"], pipe["roughness"], 10.0
    resistance = _calc_pipe_resistance_dw(diameter, roughness, viscosity, speed, _DENSITY)
    loss = get(node_fr, "head_min", -Inf) - get(node_to, "head_max", Inf)
    flow_min_loss = sign(loss) * (abs(loss) * inv(pipe["length"] * resistance))^inv(2.0)
    flow_min_dir = pipe["flow_direction"] == POSITIVE ? 0.0 : -Inf
    return max(-capacity, flow_min_loss, flow_min_dir, get(pipe, "flow_min", -Inf))
end


function _calc_pipe_flow_max_hw(
    pipe::Dict{String, <:Any}, node_fr::Dict{String, <:Any},
    node_to::Dict{String, <:Any}, capacity::Float64)
    # Calculate maximum flow bound per the Hazen-Williams head loss formulation.
    resistance = _calc_pipe_resistance_hw(pipe["diameter"], pipe["roughness"])
    loss = get(node_fr, "head_max", Inf) - get(node_to, "head_min", -Inf)
    flow_max_loss = sign(loss) * (abs(loss) * inv(pipe["length"] * resistance))^inv(1.852)
    flow_max_dir = pipe["flow_direction"] == NEGATIVE ? 0.0 : Inf
    return min(capacity, flow_max_loss, flow_max_dir, get(pipe, "flow_max", Inf))
end


function _calc_pipe_flow_max_dw(
    pipe::Dict{String, <:Any}, node_fr::Dict{String, <:Any}, node_to::Dict{String, <:Any},
    viscosity::Float64, capacity::Float64)
    # Calculate maximum flow bound per the Hazen-Williams head loss formulation.
    diameter, roughness, speed = pipe["diameter"], pipe["roughness"], 10.0
    resistance = _calc_pipe_resistance_dw(diameter, roughness, viscosity, speed, _DENSITY)
    loss = get(node_fr, "head_max", Inf) - get(node_to, "head_min", -Inf)
    flow_max_loss = sign(loss) * (abs(loss) * inv(pipe["length"] * resistance))^inv(2.0)
    flow_max_dir = pipe["flow_direction"] == NEGATIVE ? 0.0 : Inf
    return min(capacity, flow_max_loss, flow_max_dir, get(pipe, "flow_max", Inf))
end


function _get_exponent_from_head_loss_form(head_loss_form::String)
    return uppercase(head_loss_form) == "H-W" ? 1.852 : 2.0
end


"""
Computes the resistance for a pipe governed by the Hazen-Williams relationship.
"""
function _calc_pipe_resistance_hw(diameter::Float64, roughness::Float64)
    return 7.8828 * inv(0.849^1.852 * roughness^1.852 * diameter^4.8704)
end


"""
Computes the resistance for a pipe governed by the Darcy-Weisbach relationship.
"""
function _calc_pipe_resistance_dw(diameter::Float64, roughness::Float64, viscosity::Float64, speed::Float64, density::Float64)
    # Compute the Reynold's number of the fluid.
    reynolds_number = density * speed * diameter * inv(viscosity)

    # Use the same Colebrook approximation as in EPANET.
    w = 0.25 * pi * reynolds_number
    y_1 = 4.61841319859 * inv(w^0.9)
    y_2 = (roughness * inv(diameter)) * inv(3.7 * diameter) + y_1
    y_3 = -8.685889638e-01 * log(y_2)
    return 0.0826 * inv(diameter^5) * inv(y_3*y_3)
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


function is_des_pipe(pipe::Pair{Int64, <:Any})
    return any([x in ["diameters", "resistances"] for x in keys(pipe.second)])
end


function is_des_pipe(pipe::Pair{String, <:Any})
    return any([x in ["diameters", "resistances"] for x in keys(pipe.second)])
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
                r = _calc_pipe_resistance_hw(entry["diameter"], pipe["roughness"])
                resistances[a] = vcat(resistances[a], r)
            end

            resistances[a] = sort(resistances[a], rev=true)
        else
            r = _calc_pipe_resistance_hw(pipe["diameter"], pipe["roughness"])
            resistances[a] = vcat(resistances[a], r)
        end
    end

    return resistances
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
                r = _calc_pipe_resistance_dw(diameter, roughness, viscosity, 10.0, _DENSITY)
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
            r = _calc_pipe_resistance_dw(diameter, roughness, viscosity, 10.0, _DENSITY)
            resistances[a] = vcat(resistances[a], r)
        end
    end

    return resistances
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


function calc_resistance_costs_hw(pipes::Dict{Int, <:Any})
    # Create placeholder costs dictionary.
    costs = Dict([(a, Array{Float64, 1}()) for a in keys(pipes)])

    for (a, pipe) in pipes
        if haskey(pipe, "diameters")
            resistances = Array{Float64, 1}()

            for entry in pipe["diameters"]
                resistance = _calc_pipe_resistance_hw(entry["diameter"], pipe["roughness"])
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
                resistance = _calc_pipe_resistance_dw(diameter, roughness, viscosity, 10.0, _DENSITY)
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
