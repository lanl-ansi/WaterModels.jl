########################################################################
# This file defines commonly-used constraints for water systems models.
########################################################################

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

function constraint_pump_control_tank(wm::AbstractWaterModel, n_1::Int, n_2::Int, a::Int, i::Int, lt::Float64, gt::Float64, elevation::Float64)
    h = var(wm, n_2, :h, i)
    x_p = var(wm, n_2, :x_pump, a)
    x_lt = var(wm, n_2, :x_thrs_lt, a)
    x_gt = var(wm, n_2, :x_thrs_gt, a)
    x_bt = var(wm, n_2, :x_thrs_bt, a)
    x_p_tm1 = var(wm, n_1, :x_pump, a)

    h_ub = JuMP.upper_bound(h) - elevation
    h_lb = JuMP.lower_bound(h) - elevation

    # Pump constraints.
    # TODO: Clean everything below. There are probably plenty of simplifications.
    c_1 = JuMP.@constraint(wm.model, h - elevation <= x_p * gt + (1.0 - x_p) * h_ub)
    c_2 = JuMP.@constraint(wm.model, h - elevation >= (1.0 - x_p) * lt + x_p * h_lb)
    c_3 = JuMP.@constraint(wm.model, h - elevation <= x_lt * lt + (1.0 - x_lt) * h_ub)
    c_4 = JuMP.@constraint(wm.model, h - elevation >= x_gt * gt + (1.0 - x_gt) * h_lb)
    c_5 = JuMP.@constraint(wm.model, h - elevation <= x_bt * gt + (1.0 - x_bt) * h_ub)
    c_6 = JuMP.@constraint(wm.model, h - elevation >= x_bt * lt + (1.0 - x_bt) * h_lb)
    c_7 = JuMP.@constraint(wm.model, x_lt + x_gt + x_bt == 1)

    w_1 = JuMP.@variable(wm.model, binary=true, start=1)
    c_8 = JuMP.@constraint(wm.model, w_1 <= x_p_tm1)
    c_9 = JuMP.@constraint(wm.model, w_1 <= x_bt)
    c_10 = JuMP.@constraint(wm.model, w_1 >= x_bt + x_p_tm1 - 1.0)

    w_2 = JuMP.@variable(wm.model, binary=true, start=0)
    c_11 = JuMP.@constraint(wm.model, w_2 <= x_lt + x_gt)
    c_12 = JuMP.@constraint(wm.model, w_2 <= x_p)
    c_13 = JuMP.@constraint(wm.model, w_2 >= (x_lt + x_gt) + x_p - 1.0)
    c_14 = JuMP.@constraint(wm.model, w_1 + w_2 == x_p)
end

function constraint_pump_control_initial(wm::AbstractWaterModel, n::Int, a::Int, status::Bool)
    x_pump = var(wm, n, :x_pump, a)
    c = JuMP.@constraint(wm.model, x_pump == status)
end

function constraint_tank_state_initial(wm::AbstractWaterModel, n::Int, i::Int, V_0::Float64, time_step::Float64)
    V = var(wm, n, :V, i)
    c = JuMP.@constraint(wm.model, V == V_0)
    con(wm, n, :tank_state)[i] = c
end

function constraint_tank_state(wm::AbstractWaterModel, n_1::Int, n_2::Int, i::Int, time_step::Float64)
    qt = var(wm, n_1, :qt, i) # Tank outflow.
    V_1, V_2 = [var(wm, n_1, :V, i), var(wm, n_2, :V, i)]
    c = JuMP.@constraint(wm.model, V_1 - time_step * qt == V_2)
    con(wm, n_2, :tank_state)[i] = c
end

function constraint_recover_volume(wm::AbstractWaterModel, i::Int, n_1::Int, n_f::Int)
    _initialize_con_dict(wm, :recover_volume, nw=n_f)
    V_1, V_f = [var(wm, n_1, :V, i), var(wm, n_f, :V, i)]
    c = JuMP.@constraint(wm.model, V_f >= V_1)
    con(wm, n_f, :recover_volume)[i] = c
end
