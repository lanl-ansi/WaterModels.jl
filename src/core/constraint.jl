########################################################################
# This file defines commonly-used constraints for water systems models.
########################################################################

function constraint_source_head(wm::AbstractWaterModel, n::Int, i::Int, h_src::Float64)
    h = var(wm, n, :h, i)
    con(wm, n, :source_head)[i] = JuMP.@constraint(wm.model, h == h_src)
end

function constraint_flow_conservation(wm::AbstractWaterModel, n::Int, i::Int,
    a_fr::Array{Tuple{Int,Int,Int}}, a_to::Array{Tuple{Int,Int,Int}},
    reservoirs::Array{Int}, tanks::Array{Int}, demands::Dict{Int,Float64})
    q = var(wm, n, :q)
    qr = var(wm, n, :qr)
    qt = var(wm, n, :qt)

    c = JuMP.@constraint(wm.model, sum(q[a] for (a, f, t) in a_to) -
        sum(q[a] for (a, f, t) in a_fr) == -sum(qr[id] for id in reservoirs) -
        sum(qt[id] for id in tanks) + sum(demand for (id, demand) in demands))

    con(wm, n, :flow_conservation)[i] = c
end

function constraint_link_volume(wm::AbstractWaterModel, n::Int, i::Int, elevation::Float64, surface_area::Float64)
    h = var(wm, n, :h, i)
    V = var(wm, n, :V, i)
    c = JuMP.@constraint(wm.model, h - elevation == V * inv(surface_area))
    con(wm, n, :link_volume)[i] = c
end

#function constraint_head_loss_ub_pipe(wm::AbstractWaterModel, n::Int, a::Int, alpha, len, r_max)
#    L = len
#    r = r_max
#
#    dhp = var(wm, n, :dhp, a)
#    qp = var(wm, n, :qp, a)
#
#    qp_ub = JuMP.upper_bound(qp)
#    rhs_p = r * qp_ub^(alpha - 1.0) * qp
#    con_p = JuMP.@constraint(wm.model, inv(L) * dhp - rhs_p <= 0.0)
#    con(wm, n, :directed_head_loss_ub_pipe_p)[a] = con_p
#
#    dhn = var(wm, n, :dhn, a)
#    qn = var(wm, n, :qn, a)
#
#    qn_ub = JuMP.upper_bound(qn)
#    rhs_n = r * qn_ub^(alpha - 1.0) * qn
#    con_n = JuMP.@constraint(wm.model, inv(L) * dhn - rhs_n <= 0.0)
#    con(wm, n, :directed_head_loss_ub_pipe_n)[a] = con_n
#end

#function constraint_link_flow_ne(wm::AbstractUndirectedFlowModel, n::Int, a::Int)
#    if !haskey(con(wm, n), :link_undirected_flow_ne)
#        con(wm, n)[:link_undirected_flow_ne] = Dict{Int, JuMP.ConstraintRef}()
#    end
#
#    qne = var(wm, n, :qne, a)
#    q = var(wm, n, :q, a)
#
#    c = JuMP.@constraint(wm.model, sum(qne) == q)
#    con(wm, n, :link_undirected_flow_ne)[a] = c
#end
#
#function constraint_link_flow_ne(wm::AbstractDirectedFlowModel, n::Int, a::Int)
#    if !haskey(con(wm, n), :link_directed_flow_n_ne)
#        con(wm, n)[:link_directed_flow_p_ne] = Dict{Int, JuMP.ConstraintRef}()
#        con(wm, n)[:link_directed_flow_n_ne] = Dict{Int, JuMP.ConstraintRef}()
#    end
#
#    qp = var(wm, n, :qp, a)
#    qn = var(wm, n, :qn, a)
#
#    qp_ne = var(wm, n, :qp_ne, a)
#    qn_ne = var(wm, n, :qn_ne, a)
#
#    con_p = JuMP.@constraint(wm.model, sum(qp_ne) == qp)
#    con(wm, n, :link_directed_flow_p_ne)[a] = con_p
#
#    con_n = JuMP.@constraint(wm.model, sum(qn_ne) == qn)
#    con(wm, n, :link_directed_flow_n_ne)[a] = con_n
#end

function constraint_check_valve(wm::AbstractWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int)
    # Collect variables.
    q = var(wm, n, :q, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)
    x_cv = var(wm, n, :x_cv, a)

    # If the check valve is open, flow must be appreciably nonnegative.
    c_1 = JuMP.@constraint(wm.model, q <= JuMP.upper_bound(q) * x_cv)
    c_2 = JuMP.@constraint(wm.model, q >= 6.31465679e-6 * x_cv)

    ## TODO: These constraints seem to result in infeasibility in multiperiod Richmond case.
    #dh_lb = JuMP.lower_bound(h_i) - JuMP.upper_bound(h_j)
    #dh_ub = JuMP.upper_bound(h_i) - JuMP.lower_bound(h_j)
    #c_3 = JuMP.@constraint(wm.model, h_i - h_j >= (1.0 - x_cv) * dh_lb)
    #c_4 = JuMP.@constraint(wm.model, h_i - h_j <= 0.0)

    append!(con(wm, n, :check_valve)[a], [c_1, c_2, c_3, c_4])
end

# do nothing by default
function constraint_sink_flow(wm::AbstractWaterModel, n::Int, i::Int, links) end

# do nothing by default
function constraint_source_flow(wm::AbstractWaterModel, n::Int, i::Int, links) end

# no generic implementation available
function constraint_head_gain_pump(wm::AbstractWaterModel, n::Int, a::Int) end

# no generic implementation available
#function constraint_head_gain_pump_quadratic_fit(wm::AbstractWaterModel, n::Int, i::Int)
#end

function constraint_pump_control_tank(wm::AbstractWaterModel, n_1::Int,
    n_2::Int, a::Int, i::Int, lt::Float64, gt::Float64, elevation::Float64)
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
    qt = var(wm, n_1, :qt, i)
    V_1 = var(wm, n_1, :V, i)
    V_2 = var(wm, n_2, :V, i)

    c = JuMP.@constraint(wm.model, V_2 - V_1 + time_step * qt == 0.0)
    con(wm, n_2, :tank_state)[i] = c
end
