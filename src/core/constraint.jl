#########################################################################
# This file defines commonly-used constraints for water systems models. #
#########################################################################

function constraint_source_head(wm::AbstractWaterModel, n::Int, i::Int, h_s::Float64)
    c = JuMP.@constraint(wm.model, var(wm, n, :h, i) == h_s)
    con(wm, n, :source_head)[i] = c
end

function constraint_flow_conservation(wm::AbstractWaterModel, n::Int, i::Int, a_fr::Array{Tuple{Int,Int,Int}}, a_to::Array{Tuple{Int,Int,Int}}, reservoirs::Array{Int}, tanks::Array{Int}, demands::Dict{Int,Float64})
    q, qr, qt = [var(wm, n, :q), var(wm, n, :qr), var(wm, n, :qt)]

    c = JuMP.@constraint(wm.model, sum(q[a] for (a, f, t) in a_to)
        - sum(q[a] for (a, f, t) in a_fr) == -sum(qr[id] for id in reservoirs)
        - sum(qt[id] for id in tanks) + sum(demand for (id, demand) in demands))

    con(wm, n, :flow_conservation)[i] = c
end

function constraint_link_volume(wm::AbstractWaterModel, n::Int, i::Int, elevation::Float64, surface_area::Float64)
    h, V = [var(wm, n, :h, i), var(wm, n, :V, i)]
    c = JuMP.@constraint(wm.model, h - elevation == V * inv(surface_area))
    con(wm, n, :link_volume)[i] = c
end

function constraint_pump_control_initial(wm::AbstractWaterModel, n::Int, a::Int, status::Bool)
    z = var(wm, n, :z_pump, a)
    c = JuMP.@constraint(wm.model, z == status)
end

function constraint_tank_state_initial(wm::AbstractWaterModel, n::Int, i::Int, V_0::Float64, time_step::Float64)
    V = var(wm, n, :V, i)
    c = JuMP.@constraint(wm.model, V == V_0)
    con(wm, n, :tank_state)[i] = c
end

function constraint_tank_state(wm::AbstractWaterModel, n_1::Int, n_2::Int, i::Int, time_step::Float64)
    qt = var(wm, n_1, :qt, i) # Tank outflow.
    V_1, V_2 = var(wm, n_1, :V, i), var(wm, n_2, :V, i)
    c = JuMP.@constraint(wm.model, V_1 - time_step * qt == V_2)
    con(wm, n_2, :tank_state)[i] = c
end

function constraint_recover_volume(wm::AbstractWaterModel, i::Int, n_1::Int, n_f::Int)
    _initialize_con_dict(wm, :recover_volume, nw=n_f)
    V_1, V_f = var(wm, n_1, :V, i), var(wm, n_f, :V, i)
    c = JuMP.@constraint(wm.model, V_f >= V_1)
    con(wm, n_f, :recover_volume)[i] = c
end
