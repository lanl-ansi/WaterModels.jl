"""
Computes the resistance for a pipe governed by the Hazen-Williams relationship.
"""
function _calc_pipe_resistance_hazen_williams(diameter::Float64, roughness::Float64)
    return 7.8828 * inv(0.849^1.852 * roughness^1.852 * diameter^4.8704)
end


"""
Computes the resistance for a pipe governed by the Darcy-Weisbach relationship.
"""
function _calc_pipe_resistance_darcy_weisbach(diameter::Float64, roughness::Float64, viscosity::Float64, speed::Float64, density::Float64)
    # Compute the Reynold's number of the fluid.
    reynolds_number = density * speed * diameter * inv(viscosity)

    # Use the same Colebrook approximation as in EPANET.
    w = 0.25 * pi * reynolds_number
    y1 = 4.61841319859 * inv(w^0.9)
    y2 = (roughness * inv(diameter)) * inv(3.7 * diameter) + y1
    y3 = -8.685889638e-01 * log(y2)
    return 0.0826 * inv(diameter^5) * inv(y3*y3)
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
                r = _calc_pipe_resistance_hazen_williams(entry["diameter"], pipe["roughness"])
                resistances[a] = vcat(resistances[a], r)
            end

            resistances[a] = sort(resistances[a], rev=true)
        else
            r = _calc_pipe_resistance_hazen_williams(pipe["diameter"], pipe["roughness"])
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
                r = _calc_pipe_resistance_darcy_weisbach(diameter, roughness, viscosity, 10.0, 1000.0)
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
            r = _calc_pipe_resistance_darcy_weisbach(diameter, roughness, viscosity, 10.0, 1000.0)
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
                resistance = _calc_pipe_resistance_hazen_williams(entry["diameter"], pipe["roughness"])
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
                resistance = _calc_pipe_resistance_darcy_weisbach(diameter, roughness, viscosity, 10.0, 1000.0)
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
