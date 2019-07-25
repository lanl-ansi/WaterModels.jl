########################################################################
# This file defines commonly-used constraints for water systems models.
########################################################################

function constraint_flow_conservation(wm::GenericWaterModel{T}, n::Int, i::Int, node_arcs_fr, node_arcs_to, node_reservoirs, node_demands) where T <: AbstractUndirectedFlowFormulation
    q = var(wm, n, :q)
    q_r = var(wm, n, :q_r)
    q_t = var(wm, n, :q_t)

    # TBD
    # for tid in ref(wm, n, :node_tanks, i)
    #     tank = ref(wm, n, :tanks, tid)
    #     # TODO add tank vars as loads
    # end

    # Add the flow conservation constraint.
    con(wm, n, :flow_conservation)[i] = JuMP.@constraint(wm.model,
        sum(-q[l] for (l,f,t) in node_arcs_fr) +
        sum(q[l] for (l,f,t) in node_arcs_to)
        ==
        sum(-q_r[rid] for rid in node_reservoirs) +
        sum(q_t[tid] for tid in ref(wm, n, :node_tanks, i)) +
        sum(demand for (jid, demand) in node_demands)
    )
end


function constraint_flow_conservation(wm::GenericWaterModel{T}, n::Int, i::Int, node_arcs_fr, node_arcs_to, node_reservoirs, node_demands) where T <: AbstractDirectedFlowFormulation
    qn = var(wm, n, :qn)
    qp = var(wm, n, :qp)
    q_r = var(wm, n, :q_r)
    q_t = var(wm, n, :q_t)

    # TBD
    # for tid in ref(wm, n, :node_tanks, i)
    #     tank = ref(wm, n, :tanks, tid)
    #     # TODO add tank vars as loads
    # end

    # Add the flow conservation constraint.
    con(wm, n, :flow_conservation)[i] = JuMP.@constraint(wm.model,
        sum(  qp[l] - qn[l] for (l,f,t) in node_arcs_to) +
        sum( -qp[l] + qn[l] for (l,f,t) in node_arcs_fr)
        ==
        sum(-q_r[rid] for rid in node_reservoirs) +
        sum(q_t[tid] for tid in ref(wm, n, :node_tanks, i)) +
        sum(demand for (jid, demand) in node_demands)
    )
end


function constraint_resistance_selection_ne(wm::GenericWaterModel{T}, n::Int, a::Int) where T <: AbstractDirectedFlowFormulation
    if !haskey(con(wm, n), :resistance_selection_sum)
        con(wm, n)[:resistance_selection_sum] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:resistance_selection_p] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        con(wm, n)[:resistance_selection_n] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    con_sum = JuMP.@constraint(wm.model, sum(var(wm, n, :x_res, a)) == 1.0)
    con(wm, n, :resistance_selection_sum)[a] = con_sum

    con(wm, n, :resistance_selection_p)[a] = Dict{Int, JuMP.ConstraintRef}()
    con(wm, n, :resistance_selection_n)[a] = Dict{Int, JuMP.ConstraintRef}()

    for r in 1:length(ref(wm, n, :resistance, a))
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

function constraint_resistance_selection_ne(wm::GenericWaterModel{T}, n::Int, a::Int) where T <: AbstractUndirectedFlowFormulation
    if !haskey(con(wm, n), :resistance_selection_sum)
        con(wm, n)[:resistance_selection_sum] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:resistance_selection_lb] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        con(wm, n)[:resistance_selection_ub] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    con_sum = JuMP.@constraint(wm.model, sum(var(wm, n, :x_res, a)) == 1.0)
    con(wm, n, :resistance_selection_sum)[a] = con_sum

    con(wm, n, :resistance_selection_lb)[a] = Dict{Int, JuMP.ConstraintRef}()
    con(wm, n, :resistance_selection_ub)[a] = Dict{Int, JuMP.ConstraintRef}()

    for r in 1:length(ref(wm, n, :resistance, a))
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
    if !haskey(con(wm, n), :flow_direction_selection_n)
        con(wm, n)[:flow_direction_selection_n] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:flow_direction_selection_p] = Dict{Int, JuMP.ConstraintRef}()
    end

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
function constraint_flow_direction_selection_ne(wm::GenericWaterModel, n::Int, a::Int)
end

function constraint_flow_direction_selection_ne(wm::GenericWaterModel{T}, n::Int, a::Int) where T <: AbstractDirectedFlowFormulation
    if !haskey(con(wm, n), :flow_direction_selection_ne_n)
        con(wm, n)[:flow_direction_selection_ne_p] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        con(wm, n)[:flow_direction_selection_ne_n] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    con(wm, n, :flow_direction_selection_ne_p)[a] = Dict{Int, JuMP.ConstraintRef}()
    con(wm, n, :flow_direction_selection_ne_n)[a] = Dict{Int, JuMP.ConstraintRef}()

    for r in 1:length(ref(wm, n, :resistance, a))
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


function constraint_head_difference(wm::GenericWaterModel, n::Int, a::Int)
    if !haskey(con(wm, n), :head_difference_1)
        con(wm, n)[:head_difference_1] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:head_difference_2] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:head_difference_3] = Dict{Int, JuMP.ConstraintRef}()
    end

    i = ref(wm, n, :links, a)["f_id"]

    if i in collect(ids(wm, n, :reservoirs))
        h_i = ref(wm, n, :reservoirs, i)["head"]
    else
        h_i = var(wm, n, :h, i)
    end

    j = ref(wm, n, :links, a)["t_id"]

    if j in collect(ids(wm, n, :reservoirs))
        h_j = ref(wm, n, :reservoirs, j)["head"]
    else
        h_j = var(wm, n, :h, j)
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

    con_3 = JuMP.@constraint(wm.model, (h_i - h_j) - (dhp - dhn) == 0.0)
    con(wm, n, :head_difference_3)[a] = con_3
end

function constraint_potential_loss_ub_pipe_ne(wm::GenericWaterModel, n::Int, a::Int) where T <: AbstractDirectedFlowFormulation
    if !haskey(con(wm, n), :directed_potential_loss_ub_pipe_ne_n)
        con(wm, n)[:directed_potential_loss_ub_pipe_ne_n] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:directed_potential_loss_ub_pipe_ne_p] = Dict{Int, JuMP.ConstraintRef}()
    end

    alpha = ref(wm, n, :alpha)
    L = ref(wm, n, :pipes, a)["length"]
    resistances = ref(wm, n, :resistance, a)

    dhp = var(wm, n, :dhp, a)
    qp_ne_ub = JuMP.upper_bound.(var(wm, n, :qp_ne, a))
    slopes_p = resistances .* qp_ne_ub.^(alpha - 1.0)
    rhs_p = sum(slopes_p .* var(wm, n, :qp_ne, a))
    con_p = JuMP.@constraint(wm.model, inv(L) * dhp - rhs_p <= 0.0)
    con(wm, n, :directed_potential_loss_ub_pipe_ne_p)[a] = con_p

    dhn = var(wm, n, :dhn, a)
    qn_ne_ub = JuMP.upper_bound.(var(wm, n, :qn_ne, a))
    slopes_n = resistances .* qn_ne_ub.^(alpha - 1.0)
    rhs_n = sum(slopes_n .* var(wm, n, :qn_ne, a))
    con_n = JuMP.@constraint(wm.model, inv(L) * dhn - rhs_n <= 0.0)
    con(wm, n, :directed_potential_loss_ub_pipe_ne_n)[a] = con_n
end

function constraint_potential_loss_ub_pipe(wm::GenericWaterModel, n::Int, a::Int)
    if !haskey(con(wm, n), :directed_potential_loss_ub_pipe_n)
        con(wm, n)[:directed_potential_loss_ub_pipe_p] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:directed_potential_loss_ub_pipe_n] = Dict{Int, JuMP.ConstraintRef}()
    end

    alpha = ref(wm, n, :alpha)
    L = ref(wm, n, :pipes, a)["length"]
    r = maximum(ref(wm, n, :resistance, a))

    dhp = var(wm, n, :dhp, a)
    qp_ub = JuMP.upper_bound(var(wm, n, :qp, a))
    rhs_p = r * qp_ub^(alpha - 1.0) * var(wm, n, :qp, a)
    con_p = JuMP.@constraint(wm.model, inv(L) * dhp - rhs_p <= 0.0)
    con(wm, n, :directed_potential_loss_ub_pipe_p)[a] = con_p

    dhn = var(wm, n, :dhn, a)
    qn_ub = JuMP.upper_bound(var(wm, n, :qn, a))
    rhs_n = r * qn_ub^(alpha - 1.0) * var(wm, n, :qn, a)
    con_n = JuMP.@constraint(wm.model, inv(L) * dhn - rhs_n <= 0.0)
    con(wm, n, :directed_potential_loss_ub_pipe_n)[a] = con_n
end

function constraint_link_volume(wm::GenericWaterModel, i::Int, n::Int=wm.cnw)
    if !haskey(con(wm, n), :link_volume)
        con(wm, n)[:link_volume] = Dict{Int, JuMP.ConstraintRef}()
    end

    h = var(wm, n, :h, i)
    V = var(wm, n, :V, i)
    elevation = ref(wm, n, :nodes, i)["elevation"]
    surface_area = 0.25 * pi * ref(wm, n, :tanks, i)["diameter"]^2

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


# do nothing by default
function constraint_sink_flow(wm::GenericWaterModel, n::Int, i::Int)
end

"Constraint to ensure at least one direction is set to take flow to a junction with demand."
function constraint_sink_flow(wm::GenericWaterModel{T}, n::Int, i::Int) where T <: AbstractDirectedFlowFormulation
    # Collect the required variables.
    x_dir = var(wm, n, :x_dir)
    out_arcs = filter(a -> i == a.second["f_id"], ref(wm, n, :links))
    out = Array{JuMP.VariableRef}([x_dir[a] for a in keys(out_arcs)])
    in_arcs = filter(a -> i == a.second["t_id"], ref(wm, n, :links))
    in = Array{JuMP.VariableRef}([x_dir[a] for a in keys(in_arcs)])

    # Add the sink flow direction constraint.
    c = JuMP.@constraint(wm.model, sum(in) - sum(out) >= 1.0 - length(out))
    con(wm, n, :directed_sink_flow)[i] = c
end


function constraint_check_valve(wm::GenericWaterModel, a::Int, n::Int=wm.cnw)
    if !haskey(con(wm, n), :check_valve_1)
        con(wm, n)[:check_valve_1] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:check_valve_2] = Dict{Int, JuMP.ConstraintRef}()
    end

    q = var(wm, n, :q, a)
    h_i = var(wm, n, :h, ref(wm, n, :links, a)["f_id"])
    h_j = var(wm, n, :h, ref(wm, n, :links, a)["t_id"])
    x_cv = var(wm, n, :x_cv, a)

    con(wm, n, :check_valve_1)[a] = JuMP.@constraint(wm.model, q <= JuMP.upper_bound(q) * x_cv)
    #con(wm, n, :check_valve_2)[a] = JuMP.@NLconstraint(wm.model, (1 - x_cv) * h_i <= (1 - x_cv) * h_j)
end

# do nothing by default
function constraint_source_flow(wm::GenericWaterModel, n::Int, i::Int)
end

"Constraint to ensure at least one direction is set to take flow away from a source."
function constraint_source_flow(wm::GenericWaterModel{T}, n::Int, i::Int) where T <: AbstractDirectedFlowFormulation
    # Collect the required variables.
    x_dir = var(wm, n, :x_dir)
    out_arcs = filter(a -> i == a.second["f_id"], ref(wm, n, :links))
    out = Array{JuMP.VariableRef}([x_dir[a] for a in keys(out_arcs)])
    in_arcs = filter(a -> i == a.second["t_id"], ref(wm, n, :links))
    in = Array{JuMP.VariableRef}([x_dir[a] for a in keys(in_arcs)])

    # Add the source flow direction constraint.
    c = JuMP.@constraint(wm.model, sum(out) - sum(in) >= 1.0 - length(in))
    con(wm, n, :directed_source_flow)[i] = c
end



# do nothing by default
function constraint_potential_loss_pump(wm::GenericWaterModel, n::Int, i::Int)
end

# do nothing by default
function constraint_head_gain_pump_quadratic_fit(wm::GenericWaterModel, n::Int, i::Int)
end

""
function constraint_tank_state_initial(wm::GenericWaterModel, n::Int, i::Int, initial_volume::Float64, time_step::Float64)
    if !haskey(con(wm, n), :tank_state_1)
        con(wm, n)[:tank_state_1] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:tank_state_2] = Dict{Int, JuMP.ConstraintRef}()
    end

    V = var(wm, n, :V, i)
    q_t = var(wm, n, :q_t, i)

    #con(wm, n, :tank_state_1)[i] = JuMP.@constraint(wm.model, V - q_t * time_step == initial_volume)
    #con(wm, n, :tank_state_1)[i] = JuMP.@constraint(wm.model, V == initial_volume)
    #con(wm, n, :tank_state_2)[i] = JuMP.@constraint(wm.model, q_t == 0.0)
end

function constraint_tank_state(wm::GenericWaterModel, n_1::Int, n_2::Int, i::Int, time_step::Float64)
    if !haskey(con(wm, n), :tank_state_1)
        con(wm, n)[:tank_state_1] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:tank_state_2] = Dict{Int, JuMP.ConstraintRef}()
    end

    V_1 = var(wm, n_1, :V, i)
    V_2 = var(wm, n_2, :V, i)
    q_t = var(wm, n_2, :q_t, i)

    con(wm, n, :tank_state_1)[i] = JuMP.@constraint(wm.model, V_2 - V_1 == time_step * q_t)
end
