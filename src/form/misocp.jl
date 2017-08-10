# Define MISOCP implementations of Water Models

export
    MISOCPWaterModel, StandardMISOCPForm

abstract AbstractMISOCPForm <: AbstractWaterFormulation

type StandardMISOCPForm <: AbstractMISOCPForm end
typealias MISOCPWaterModel GenericWaterModel{StandardMISOCPForm}

# default MISOCP constructor
function MISOCPWaterModel(data::Dict{AbstractString,Any}; kwargs...)
    return GenericWaterModel(data, StandardMISOCPForm(); kwargs...)
end

# variables associated with the flux squared
function variable_flux_square{T <: AbstractMISOCPForm}(gm::GenericWaterModel{T})
    max_flow = wm.data["max_flow"]
    @variable(wm.model, 0 <= l[i in [wm.set.pipe_indexes; wm.set.resistor_indexes]] <= 1/wm.set.connections[i]["resistance"] * max_flow^2, start = getstart(wm.set.connections, i, "l_start", 0))
    return l
end

# # variables associated with the flux squared
# function variable_flux_square_ne{T <: AbstractMISOCPForm}(wm::GenericWaterModel{T})
#     max_flow = wm.data["max_flow"]
#     @variable(wm.model, 0 <= l_ne[i in wm.set.new_pipe_indexes] <= 1/wm.set.new_connections[i]["resistance"] * max_flow^2, start = getstart(wm.set.new_connections, i, "l_start", 0))
#     return l_ne
# end

#Weymouth equation with discrete direction variables
function constraint_weisbach{T <: AbstractMISOCPForm}(wm::GenericWaterModel{T}, pipe)

    pipe_idx = pipe["index"]
    i_junction_idx = pipe["f_junction"]
    j_junction_idx = pipe["t_junction"]

    hi = getvariable(wm.model, :h_water)[i_junction_idx]
    hj = getvariable(wm.model, :h_water)[j_junction_idx]
    yp = getvariable(wm.model, :yp)[pipe_idx]
    yn = getvariable(wm.model, :yn)[pipe_idx]
    l  = getvariable(wm.model, :l)[pipe_idx]
    f  = getvariable(wm.model, :f)[pipe_idx]

    hd_max = pipe["hd_max"] #i["pmax"]^2 - j["pmin"]^2;
    hd_min = pipe["hd_min"] # i["pmin"]^2 - j["pmax"]^2;
    max_flow = wm.data["max_flow"]

    c1 = @constraint(wm.model, l >= hj - hi + hd_min*(yp - yn + 1))
    c2 = @constraint(wm.model, l >= hi - hj + hd_max*(yp - yn - 1))
    c3 = @constraint(wm.model, l <= hj - hi + hd_max*(yp - yn + 1))
    c4 = @constraint(wm.model, l <= hi - hj + hd_min*(yp - yn - 1))
    c5 = @constraint(wm.model, pipe["resistance"]*l >= f^2)

    return Set([c1, c2, c3, c4, c5])
end

#Weymouth equation with fixed direction
function constraint_weisbach_fixed_direction{T <: AbstractMISOCPForm}(wm::GenericWaterModel{T}, pipe)

    pipe_idx = pipe["index"]
    i_junction_idx = pipe["f_junction"]
    j_junction_idx = pipe["t_junction"]

    hi = getvariable(wm.model, :h_water)[i_junction_idx]
    hj = getvariable(wm.model, :h_water)[j_junction_idx]
    yp = pipe["yp"]
    yn = pipe["yn"]
    l  = getvariable(wm.model, :l)[pipe_idx]
    f  = getvariable(wm.model, :f)[pipe_idx]

    hd_max = pipe["hd_max"] #i["pmax"]^2 - j["pmin"]^2;
    hd_min = pipe["hd_min"] # i["pmin"]^2 - j["pmax"]^2;
    max_flow = wm.data["max_flow"]

    c1 = @constraint(wm.model, l >= hj - hi + hd_min*(yp - yn + 1))
    c2 = @constraint(wm.model, l >= hi - hj + hd_max*(yp - yn - 1))
    c3 = @constraint(wm.model, l <= hj - hi + hd_max*(yp - yn + 1))
    c4 = @constraint(wm.model, l <= hi - hj + hd_min*(yp - yn - 1))
    c5 = @constraint(wm.model, pipe["resistance"]*l >= f^2)

    return Set([c1, c2, c3, c4, c5])
end
