"Create network expansion flow variables for directed flow formulations."
function variable_head(wm::MICPEWaterModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    # Initialize variables associated with heads.
    h = var(wm, nw)[:h] = JuMP.@variable(wm.model,
        [i in ids(wm, nw, :node)], base_name="$(nw)_h",
        start=comp_start_value(ref(wm, nw, :node, i), "h_start"))

    # Create dictionary for undirected design flow variables (dhp_ne and dhn_ne).
    dhp_ne = var(wm, nw)[:dhp_ne] = Dict{Int,Array{JuMP.VariableRef}}()
    dhn_ne = var(wm, nw)[:dhn_ne] = Dict{Int,Array{JuMP.VariableRef}}()

    # Initialize the variables. (The default start value of 1.0e-6 is crucial.)
    for a in ids(wm, nw, :link_ne)
        dhp_ne[a] = var(wm, nw, :dhp_ne)[a] = JuMP.@variable(wm.model,
            [r in 1:length(ref(wm, nw, :resistance, a))], lower_bound=0.0,
            base_name="$(nw)_dhp_ne[$(a)]",
            start=comp_start_value(ref(wm, nw, :link_ne, a), "dhp_ne_start", r, 1.0e-6))

        dhn_ne[a] = var(wm, nw, :dhn_ne)[a] = JuMP.@variable(wm.model,
            [r in 1:length(ref(wm, nw, :resistance, a))], lower_bound=0.0,
            base_name="$(nw)_dhn_ne[$(a)]",
            start=comp_start_value(ref(wm, nw, :link_ne, a), "dhn_ne_start", r, 1.0e-6))
    end

    if bounded # If the variables are bounded, apply the bounds.
        # Get the head bound variables.
        h_lb, h_ub = calc_head_bounds(wm, nw)

        # Set lower and upper bounds on heads.
        for (i, node) in ref(wm, nw, :node)
            JuMP.set_lower_bound(h[i], h_lb[i])
            JuMP.set_upper_bound(h[i], h_ub[i])
        end

        # Set lower and upper bounds on head differences.
        for (a, link) in ref(wm, nw, :link)
            i, j = [link["node_fr"], link["node_to"]]
            JuMP.set_upper_bound.(dhp_ne[a], max(0.0, h_ub[i] - h_lb[j]))
            JuMP.set_upper_bound.(dhn_ne[a], max(0.0, h_ub[j] - h_lb[i]))
        end
    end

    # Create expressions capturing the relationships among q, dhp_ne, and dhn_ne.
    var(wm, nw)[:dhp] = JuMP.@expression(wm.model,
        [a in ids(wm, nw, :link_ne)], sum(var(wm, nw, :dhp_ne, a)))

    var(wm, nw)[:dhn] = JuMP.@expression(wm.model,
        [a in ids(wm, nw, :link_ne)], sum(var(wm, nw, :dhn_ne, a)))

    # Report back head values as part of the solution.
    report && sol_component_value(wm, nw, :node, :h, ids(wm, nw, :node), h)
end

function constraint_head_difference(wm::MICPEWaterModel, n::Int, a::Int, node_fr::Int, node_to::Int, head_fr, head_to)
    x_dir = var(wm, n, :x_dir, a)

    for (r_id, r) in enumerate(ref(wm, n, :resistance, a))
        x_res = var(wm, n, :x_res, a)[r_id]
        dhp, dhn = [var(wm, n, :dhp_ne, a)[r_id], var(wm, n, :dhn_ne, a)[r_id]]
        dhp_ub, dhn_ub = [JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)]

        c_p_r = JuMP.@constraint(wm.model, dhp <= dhp_ub * x_res)
        c_p_d = JuMP.@constraint(wm.model, dhp <= dhp_ub * x_dir)
        c_n_r = JuMP.@constraint(wm.model, dhn <= dhn_ub * x_res)
        c_n_d = JuMP.@constraint(wm.model, dhn <= dhn_ub * (1.0 - x_dir))

        # Append the constraint array.
        append!(con(wm, n, :head_loss)[a], [c_p_r, c_p_d, c_n_r, c_n_d])
    end

    dhp, dhn = [var(wm, n, :dhp, a), var(wm, n, :dhn, a)]
    h_i = head_fr == nothing ? var(wm, n, :h, node_fr) : head_fr
    h_j = head_to == nothing ? var(wm, n, :h, node_to) : head_to
    c_e = JuMP.@constraint(wm.model, dhp - dhn == h_i - h_j)
    append!(con(wm, n, :head_loss)[a], [c_e])
end

function constraint_head_loss_pipe_ne(wm::MICPEWaterModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, pipe_resistances) 
    for (r_id, r) in enumerate(pipe_resistances)
        # Collect directed flow and head difference variables.
        qp, qn = [var(wm, n, :qp_ne, a)[r_id], var(wm, n, :qn_ne, a)[r_id]]
        dhp, dhn = [var(wm, n, :dhp_ne, a)[r_id], var(wm, n, :dhn_ne, a)[r_id]]

        # Add the relaxed head loss constraints.
        c_p = JuMP.@NLconstraint(wm.model, r * head_loss(qp) <= inv(L) * dhp)
        c_n = JuMP.@NLconstraint(wm.model, r * head_loss(qn) <= inv(L) * dhn)

        # Append the constraint array.
        append!(con(wm, n, :head_loss, a), [c_p, c_n])
    end
end

function constraint_energy_conservation(wm::MICPEWaterModel, n::Int, r, L, alpha)
    # Gather common variables.
    qp, qn = [var(wm, n, :qp_ne), var(wm, n, :qn_ne)]
    h, dhp, dhn = [var(wm, n, :h), var(wm, n, :dhp_ne), var(wm, n, :dhn_ne)]

    # Construct the sums used in the strong duality constraint.
    f_1 = JuMP.@NLexpression(wm.model, sum(L[a]*sum(r[a][r_id]
        * (primal_energy(qp[a][r_id]) + primal_energy(qn[a][r_id]))
        for r_id in 1:length(r[a])) for a in ids(wm, n, :pipe_ne)))

    f_2 = JuMP.@NLexpression(wm.model, sum(reservoir["head"]
        * sum(sum(qp[a][r_id] - qn[a][r_id] for r_id in
        1:length(r[a])) for (a, f, t) in ref(wm, n, :node_arc_fr, i))
        for (i, reservoir) in ref(wm, n, :reservoir)))

    f_3 = JuMP.@NLexpression(wm.model, sum(sum((L[a]*r[a][r_id])^-inv(alpha)
        * (dual_energy(dhp[a][r_id]) + dual_energy(dhn[a][r_id]))
        for r_id in 1:length(r[a])) for a in ids(wm, n, :pipe_ne)))

    f_4 = JuMP.@NLexpression(wm.model, sum(junction["demand"] * h[i]
        for (i, junction) in ref(wm, n, :junction)))

    # Add the strong duality (energy conservation) constraint.
    c = JuMP.@NLconstraint(wm.model, f_1 - f_2 + f_3 + f_4 <= 0.0)
    con(wm, n)[:energy_conservation] = c
end

function constraint_head_loss_ub_pipe_ne(wm::MICPEWaterModel, n::Int, a::Int, alpha::Float64, L::Float64, pipe_resistances)
    for (r_id, r) in enumerate(pipe_resistances)
        qp, qn = [var(wm, n, :qp_ne, a)[r_id], var(wm, n, :qn_ne, a)[r_id]]
        dhp, dhn = [var(wm, n, :dhp_ne, a)[r_id], var(wm, n, :dhn_ne, a)[r_id]]
        qp_ub, qn_ub = [JuMP.upper_bound(qp), JuMP.upper_bound(qn)]
        slope_p, slope_n = r .* [qp_ub, qn_ub].^(alpha - 1.0)

        c_p = JuMP.@constraint(wm.model, inv(L)*dhp <= slope_p * qp)
        c_n = JuMP.@constraint(wm.model, inv(L)*dhn <= slope_n * qn)
        append!(con(wm, n, :head_loss)[a], [c_p, c_n])
    end
end
