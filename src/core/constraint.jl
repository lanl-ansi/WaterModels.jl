########################################################################
# This file defines commonly-used constraints for water systems models.
########################################################################

function constraint_flow_conservation(wm::GenericWaterModel{T}, n::Int, i::Int, node_arcs_fr, node_arcs_to, node_reservoirs, node_tanks, node_demands) where T <: AbstractUndirectedFlowFormulation
    q = var(wm, n, :q)
    q_r = var(wm, n, :q_r)
    q_t = var(wm, n, :q_t)

    # Add the flow conservation constraint.
    con(wm, n, :flow_conservation)[i] = JuMP.@constraint(wm.model,
        sum(-q[l] for (l,f,t) in node_arcs_fr) +
        sum(q[l] for (l,f,t) in node_arcs_to)
        ==
        sum(-q_r[rid] for rid in node_reservoirs) +
        sum(-q_t[tid] for tid in node_tanks) +
        sum(demand for (jid, demand) in node_demands)
    )
end


function constraint_flow_conservation(wm::GenericWaterModel{T}, n::Int, i::Int, node_arcs_fr, node_arcs_to, node_reservoirs, node_tanks, node_demands) where T <: AbstractDirectedFlowFormulation
    qn = var(wm, n, :qn)
    qp = var(wm, n, :qp)
    q_r = var(wm, n, :q_r)
    q_t = var(wm, n, :q_t)

    # Add the flow conservation constraint.
    con(wm, n, :flow_conservation)[i] = JuMP.@constraint(wm.model,
        sum(  qp[l] - qn[l] for (l,f,t) in node_arcs_to) +
        sum( -qp[l] + qn[l] for (l,f,t) in node_arcs_fr)
        ==
        sum(-q_r[rid] for rid in node_reservoirs) +
        sum(-q_t[tid] for tid in node_tanks) +
        sum(demand for (jid, demand) in node_demands)
    )
end


function constraint_resistance_selection_ne(wm::GenericWaterModel{T}, n::Int, a::Int, pipe_resistances) where T <: AbstractDirectedFlowFormulation
    if !haskey(con(wm, n), :resistance_selection_sum)
        con(wm, n)[:resistance_selection_sum] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:resistance_selection_p] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        con(wm, n)[:resistance_selection_n] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    con_sum = JuMP.@constraint(wm.model, sum(var(wm, n, :x_res, a)) == 1.0)
    con(wm, n, :resistance_selection_sum)[a] = con_sum

    con(wm, n, :resistance_selection_p)[a] = Dict{Int, JuMP.ConstraintRef}()
    con(wm, n, :resistance_selection_n)[a] = Dict{Int, JuMP.ConstraintRef}()

    for r in 1:length(pipe_resistances)
        x_res = var(wm, n, :x_res, a)[r]

        qp_ne = var(wm, n, :qp_ne, a)[r]
        qp_ne_ub = JuMP.upper_bound(qp_ne)
        con_p = JuMP.@constraint(wm.model, qp_ne - qp_ne_ub * x_res <= 0.0)
        con(wm, n, :resistance_selection_p, a)[r] = con_p

        qn_ne = var(wm, n, :qn_ne, a)[r]
        qn_ne_ub = JuMP.upper_bound(qn_ne)
        con_n = JuMP.@constraint(wm.model, qn_ne - qn_ne_ub * x_res <= 0.0)
        con(wm, n, :resistance_selection_n, a)[r] = con_n
    end
end

function constraint_resistance_selection_ne(wm::GenericWaterModel{T}, n::Int, a::Int, pipe_resistances) where T <: AbstractUndirectedFlowFormulation
    if !haskey(con(wm, n), :resistance_selection_sum)
        con(wm, n)[:resistance_selection_sum] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:resistance_selection_lb] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        con(wm, n)[:resistance_selection_ub] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    con_sum = JuMP.@constraint(wm.model, sum(var(wm, n, :x_res, a)) == 1.0)
    con(wm, n, :resistance_selection_sum)[a] = con_sum

    con(wm, n, :resistance_selection_lb)[a] = Dict{Int, JuMP.ConstraintRef}()
    con(wm, n, :resistance_selection_ub)[a] = Dict{Int, JuMP.ConstraintRef}()

    for r in 1:length(pipe_resistances)
        x_res = var(wm, n, :x_res, a)[r]

        q_ne = var(wm, n, :q_ne, a)[r]
        q_ne_lb = JuMP.lower_bound(q_ne)
        con_lb = JuMP.@constraint(wm.model, q_ne - q_ne_lb * x_res >= 0.0)
        con(wm, n, :resistance_selection_lb, a)[r] = con_lb

        q_ne = var(wm, n, :q_ne, a)[r]
        q_ne_ub = JuMP.upper_bound(q_ne)
        con_ub = JuMP.@constraint(wm.model, q_ne - q_ne_ub * x_res <= 0.0)
        con(wm, n, :resistance_selection_ub, a)[r] = con_ub
    end
end


# do nothing by default
function constraint_flow_direction_selection(wm::GenericWaterModel, n::Int, a::Int)
end

function constraint_flow_direction_selection(wm::GenericWaterModel{T}, n::Int, a::Int) where T <: AbstractDirectedFlowFormulation
    x_dir = var(wm, n, :x_dir, a)

    qp = var(wm, n, :qp, a)
    qp_ub = JuMP.upper_bound(qp)
    con_p = JuMP.@constraint(wm.model, qp - qp_ub * x_dir <= 0.0)
    con(wm, n, :flow_direction_selection_p)[a] = con_p

    qn = var(wm, n, :qn, a)
    qn_ub = JuMP.upper_bound(qn)
    con_n = JuMP.@constraint(wm.model, qn - qn_ub * (1.0 - x_dir) <= 0.0)
    con(wm, n, :flow_direction_selection_n)[a] = con_n
end


# do nothing by default
function constraint_flow_direction_selection_ne(wm::GenericWaterModel, n::Int, a::Int, pipe_resistances)
end

function constraint_flow_direction_selection_ne(wm::GenericWaterModel{T}, n::Int, a::Int, pipe_resistances) where T <: AbstractDirectedFlowFormulation
    for r in 1:length(pipe_resistances)
        x_dir = var(wm, n, :x_dir, a)

        qp_ne = var(wm, n, :qp_ne, a)[r]
        qp_ne_ub = JuMP.upper_bound(qp_ne)
        con_p = JuMP.@constraint(wm.model, qp_ne - qp_ne_ub * x_dir <= 0.0)
        con(wm, n, :flow_direction_selection_ne_p, a)[r] = con_p

        qn_ne = var(wm, n, :qn_ne, a)[r]
        qn_ne_ub = JuMP.upper_bound(qn_ne)
        con_n = JuMP.@constraint(wm.model, qn_ne - qn_ne_ub * (1.0 - x_dir) <= 0.0)
        con(wm, n, :flow_direction_selection_ne_n, a)[r] = con_n
    end
end


function constraint_head_difference(wm::GenericWaterModel, n::Int, a::Int, f_id, t_id, head_fr, head_to)
    if head_fr == nothing
        head_fr = var(wm, n, :h, f_id)
    end

    if head_to == nothing
        head_to = var(wm, n, :h, t_id)
    end

    x_dir = var(wm, n, :x_dir, a)

    dhp = var(wm, n, :dhp, a)
    dhp_ub = JuMP.upper_bound(dhp)
    con_1 = JuMP.@constraint(wm.model, dhp - dhp_ub * x_dir <= 0.0)
    con(wm, n, :head_difference_1)[a] = con_1

    dhn = var(wm, n, :dhn, a)
    dhn_ub = JuMP.upper_bound(dhn)
    con_2 = JuMP.@constraint(wm.model, dhn - dhn_ub * (1.0 - x_dir) <= 0.0)
    con(wm, n, :head_difference_2)[a] = con_2

    con_3 = JuMP.@constraint(wm.model, (head_fr - head_to) - (dhp - dhn) == 0.0)
    con(wm, n, :head_difference_3)[a] = con_3
end

function constraint_potential_loss_ub_pipe_ne(wm::GenericWaterModel, n::Int, a::Int, alpha, len, pipe_resistances) where T <: AbstractDirectedFlowFormulation
    if !haskey(con(wm, n), :directed_potential_loss_ub_pipe_ne_n)
        con(wm, n)[:directed_potential_loss_ub_pipe_ne_n] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:directed_potential_loss_ub_pipe_ne_p] = Dict{Int, JuMP.ConstraintRef}()
    end

    L = len

    dhp = var(wm, n, :dhp, a)
    qp_ne = var(wm, n, :qp_ne, a)
    qp_ne_ub = JuMP.upper_bound.(qp_ne)
    slopes_p = pipe_resistances .* qp_ne_ub.^(alpha - 1.0)
    rhs_p = sum(slopes_p .* qp_ne)
    con_p = JuMP.@constraint(wm.model, inv(L) * dhp - rhs_p <= 0.0)
    con(wm, n, :directed_potential_loss_ub_pipe_ne_p)[a] = con_p

    dhn = var(wm, n, :dhn, a)
    qn_ne = var(wm, n, :qn_ne, a)
    qn_ne_ub = JuMP.upper_bound.(qn_ne)
    slopes_n = pipe_resistances .* qn_ne_ub.^(alpha - 1.0)
    rhs_n = sum(slopes_n .* qn_ne)
    con_n = JuMP.@constraint(wm.model, inv(L) * dhn - rhs_n <= 0.0)
    con(wm, n, :directed_potential_loss_ub_pipe_ne_n)[a] = con_n
end

function constraint_potential_loss_ub_pipe(wm::GenericWaterModel, n::Int, a::Int, alpha, len, r_max)
    L = len
    r = r_max

    dhp = var(wm, n, :dhp, a)
    qp = var(wm, n, :qp, a)

    qp_ub = JuMP.upper_bound(qp)
    rhs_p = r * qp_ub^(alpha - 1.0) * qp
    con_p = JuMP.@constraint(wm.model, inv(L) * dhp - rhs_p <= 0.0)
    con(wm, n, :directed_potential_loss_ub_pipe_p)[a] = con_p

    dhn = var(wm, n, :dhn, a)
    qn = var(wm, n, :qn, a)

    qn_ub = JuMP.upper_bound(qn)
    rhs_n = r * qn_ub^(alpha - 1.0) * qn
    con_n = JuMP.@constraint(wm.model, inv(L) * dhn - rhs_n <= 0.0)
    con(wm, n, :directed_potential_loss_ub_pipe_n)[a] = con_n
end


function constraint_link_volume(wm::GenericWaterModel, n::Int, i::Int, elevation, surface_area)
    h = var(wm, n, :h, i)
    V = var(wm, n, :V, i)

    c = JuMP.@constraint(wm.model, h - elevation == inv(surface_area) * V)
    con(wm, n, :link_volume)[i] = c
end


function constraint_link_flow_ne(wm::GenericWaterModel{T}, n::Int, a::Int) where T <: AbstractUndirectedFlowFormulation
    if !haskey(con(wm, n), :link_undirected_flow_ne)
        con(wm, n)[:link_undirected_flow_ne] = Dict{Int, JuMP.ConstraintRef}()
    end

    q_ne = var(wm, n, :q_ne, a)
    q = var(wm, n, :q, a)

    c = JuMP.@constraint(wm.model, sum(q_ne) == q)
    con(wm, n, :link_undirected_flow_ne)[a] = c
end

function constraint_link_flow_ne(wm::GenericWaterModel{T}, n::Int, a::Int) where T <: AbstractDirectedFlowFormulation
    if !haskey(con(wm, n), :link_directed_flow_n_ne)
        con(wm, n)[:link_directed_flow_p_ne] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:link_directed_flow_n_ne] = Dict{Int, JuMP.ConstraintRef}()
    end

    qp = var(wm, n, :qp, a)
    qn = var(wm, n, :qn, a)

    qp_ne = var(wm, n, :qp_ne, a)
    qn_ne = var(wm, n, :qn_ne, a)

    con_p = JuMP.@constraint(wm.model, sum(qp_ne) == qp)
    con(wm, n, :link_directed_flow_p_ne)[a] = con_p

    con_n = JuMP.@constraint(wm.model, sum(qn_ne) == qn)
    con(wm, n, :link_directed_flow_n_ne)[a] = con_n
end


# do nothing by default
function constraint_link_flow(wm::GenericWaterModel, n::Int, a::Int)
end

# link undirected flow variables into the directed flow variable
function constraint_link_flow(wm::GenericWaterModel{T}, n::Int, a::Int) where T <: AbstractDirectedFlowFormulation
    q = var(wm, n, :q, a)
    qp = var(wm, n, :qp, a)
    qn = var(wm, n, :qn, a)

    con(wm, n, :link_directed_flow)[a] = JuMP.@constraint(wm.model, qp - qn == q)
end


function constraint_check_valve(wm::GenericWaterModel, n::Int, a::Int, f_id::Int, t_id::Int)
    q = var(wm, n, :q, a)
    h_i = var(wm, n, :h, f_id)
    h_j = var(wm, n, :h, t_id)
    x_cv = var(wm, n, :x_cv, a)

    con(wm, n, :check_valve_1)[a] = JuMP.@constraint(wm.model, q <= JuMP.upper_bound(q) * x_cv)
    con(wm, n, :check_valve_2)[a] = JuMP.@NLconstraint(wm.model, x_cv * h_i >= x_cv * h_j)
    con(wm, n, :check_valve_3)[a] = JuMP.@NLconstraint(wm.model, (1 - x_cv) * h_i <= (1 - x_cv) * h_j)
end


# do nothing by default
function constraint_sink_flow(wm::GenericWaterModel, n::Int, i::Int, links)
end

"Constraint to ensure at least one direction is set to take flow to a junction with demand."
function constraint_sink_flow(wm::GenericWaterModel{T}, n::Int, i::Int, links) where T <: AbstractDirectedFlowFormulation
    # Collect the required variables.
    x_dir = var(wm, n, :x_dir)
    out_arcs = filter(a -> i == a.second["f_id"], links)
    out = Array{JuMP.VariableRef}([x_dir[a] for a in keys(out_arcs)])
    in_arcs = filter(a -> i == a.second["t_id"], links)
    in = Array{JuMP.VariableRef}([x_dir[a] for a in keys(in_arcs)])

    # Add the sink flow direction constraint.
    c = JuMP.@constraint(wm.model, sum(in) - sum(out) >= 1.0 - length(out))
    con(wm, n, :directed_sink_flow)[i] = c
end


# do nothing by default
function constraint_source_flow(wm::GenericWaterModel, n::Int, i::Int, links)
end

"Constraint to ensure at least one direction is set to take flow away from a source."
function constraint_source_flow(wm::GenericWaterModel{T}, n::Int, i::Int, links) where T <: AbstractDirectedFlowFormulation
    # Collect the required variables.
    x_dir = var(wm, n, :x_dir)
    out_arcs = filter(a -> i == a.second["f_id"], links)
    out = Array{JuMP.VariableRef}([x_dir[a] for a in keys(out_arcs)])
    in_arcs = filter(a -> i == a.second["t_id"], links)
    in = Array{JuMP.VariableRef}([x_dir[a] for a in keys(in_arcs)])

    # Add the source flow direction constraint.
    c = JuMP.@constraint(wm.model, sum(out) - sum(in) >= 1.0 - length(in))
    con(wm, n, :directed_source_flow)[i] = c
end

# no generic implementation available
#function constraint_potential_loss_pump(wm::GenericWaterModel, n::Int, i::Int)
#end

# no generic implementation available
#function constraint_head_gain_pump_quadratic_fit(wm::GenericWaterModel, n::Int, i::Int)
#end

function constraint_pump_control_tank(wm::GenericWaterModel, n_1::Int, n_2::Int, a::Int, i::Int, lt::Float64, gt::Float64, elevation::Float64)
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
    c_1 = JuMP.@constraint(wm.model, h - elevation <= x_p * gt + (1 - x_p) * h_ub)
    c_2 = JuMP.@constraint(wm.model, h - elevation >= (1 - x_p) * lt + x_p * h_lb)
    c_3 = JuMP.@constraint(wm.model, h - elevation <= x_lt * lt + (1 - x_lt) * h_ub)
    c_4 = JuMP.@constraint(wm.model, h - elevation >= x_gt * gt + (1 - x_gt) * h_lb)
    c_5 = JuMP.@constraint(wm.model, h - elevation <= x_bt * gt + (1 - x_bt) * h_ub)
    c_6 = JuMP.@constraint(wm.model, h - elevation >= x_bt * lt + (1 - x_bt) * h_lb)
    c_7 = JuMP.@constraint(wm.model, x_lt + x_gt + x_bt == 1)

    w_1 = JuMP.@variable(wm.model, binary=true, start=1.0)
    c_8 = JuMP.@constraint(wm.model, w_1 <= x_p_tm1)
    c_9 = JuMP.@constraint(wm.model, w_1 <= x_bt)
    c_10 = JuMP.@constraint(wm.model, w_1 >= x_bt + x_p_tm1 - 1.0)

    w_2 = JuMP.@variable(wm.model, binary=true, start=0.0)
    c_11 = JuMP.@constraint(wm.model, w_2 <= x_lt + x_gt)
    c_12 = JuMP.@constraint(wm.model, w_2 <= x_p)
    c_13 = JuMP.@constraint(wm.model, w_2 >= (x_lt + x_gt) + x_p - 1.0)
    c_14 = JuMP.@constraint(wm.model, w_1 + w_2 == x_p)
end

function constraint_pump_control_initial(wm::GenericWaterModel, n::Int, a::Int, status::Bool)
    x = var(wm, n, :x_pump, a)
    c_1 = JuMP.@constraint(wm.model, x == status)
end

""
function constraint_tank_state_initial(wm::GenericWaterModel, n::Int, i::Int, initial_volume::Float64, time_step::Float64)
    if !haskey(con(wm, n), :tank_state)
        con(wm, n)[:tank_state] = Dict{Int, JuMP.ConstraintRef}()
    end

    V = var(wm, n, :V, i)
    con(wm, n, :tank_state)[i] = JuMP.@constraint(wm.model, V == initial_volume)
end


""
function constraint_tank_state(wm::GenericWaterModel, n_1::Int, n_2::Int, i::Int, time_step::Float64)
    if !haskey(con(wm, n_2), :tank_state)
        con(wm, n_2)[:tank_state] = Dict{Int, JuMP.ConstraintRef}()
    end

    V_1 = var(wm, n_1, :V, i)
    V_2 = var(wm, n_2, :V, i)
    q_t = var(wm, n_1, :q_t, i)

    con(wm, n_2, :tank_state)[i] = JuMP.@constraint(wm.model, V_2 - V_1 + time_step * q_t == 0.0)
end
