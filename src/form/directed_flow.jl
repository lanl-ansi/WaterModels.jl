# Constraints and variables common to all formulations with directed flows. In these
# formulations, the variables qp correspond to flow from i to j, and the variables qn
# correspond to flow from j to i. That is, when qp is nonzero, qn should be zero, and when
# qn is nonzero, qp should be zero.


function _collect_directions(wm::AbstractWaterModel, n::Int, a::Int)
    comps = ["check_valve", "pipe", "pressure_reducing_valve", "pump", "shutoff_valve"]
    comps = filter(x -> a in ids(wm, n, Symbol(x)), comps)
    return sum(var(wm, n, Symbol("y_" * c), a) for c in comps)
end


"Initialize variables associated with flow direction. If this variable is equal
to one, the flow direction is from i to j. If it is equal to zero, the flow
direction is from j to i."
function _variable_component_direction(
    wm::AbstractDirectedModel, component_name::String; nw::Int=wm.cnw, report::Bool=true)
    # Store the corresponding component symbol.
    comp_sym = Symbol(component_name)

    # Initialize variables associated with positive flows.
    y = var(wm, nw)[Symbol("y_" * component_name)] = JuMP.@variable(
        wm.model, [a in ids(wm, nw, comp_sym)], binary=true, base_name="$(nw)_y",
        start=comp_start_value(ref(wm, nw, comp_sym, a), "y_start"))

    # Report back flow values as part of the solution.
    report && sol_component_value(wm, nw, comp_sym, :y, ids(wm, nw, comp_sym), y)
end


"Create head differences variables common to all directed flow models for a component."
function _variable_component_head_difference(
    wm::AbstractDirectedModel, component_name::String; nw::Int=wm.cnw, bounded::Bool=true,
    report::Bool=true)
    # Store the corresponding component symbol.
    comp_sym = Symbol(component_name)

    # Initialize directed variables associated with positive flows.
    dhp = var(wm, nw)[Symbol("dhp_" * component_name)] = JuMP.@variable(
        wm.model, [a in ids(wm, nw, comp_sym)], lower_bound=0.0, base_name="$(nw)_dhp",
        start=comp_start_value(ref(wm, nw, comp_sym, a), "dhp_start"))

    # Initialize directed variables associated with negative flows.
    dhn = var(wm, nw)[Symbol("dhn_" * component_name)] = JuMP.@variable(
        wm.model, [a in ids(wm, nw, comp_sym)], lower_bound=0.0, base_name="$(nw)_dhn",
        start=comp_start_value(ref(wm, nw, comp_sym, a), "dhn_start"))

    if bounded # Bound flow-related variables if desired.
        # Get the head bound variables.
        h_lb, h_ub = calc_head_bounds(wm, nw)

        # Set lower and upper bounds on head differences.
        for (a, comp) in ref(wm, nw, comp_sym)
            i, j = comp["node_fr"], comp["node_to"]
            JuMP.set_upper_bound(dhp[a], max(0.0, h_ub[i] - h_lb[j]))
            JuMP.set_upper_bound(dhn[a], max(0.0, h_ub[j] - h_lb[i]))
        end
    end

    # Report positive head difference values as part of the solution.
    report && sol_component_value(wm, nw, comp_sym, :dhp, ids(wm, nw, comp_sym), dhp)

    # Report positive head difference values as part of the solution.
    report && sol_component_value(wm, nw, comp_sym, :dhn, ids(wm, nw, comp_sym), dhn)
end


"Create flow variables that are common to all directed flow models for a component."
function _variable_component_flow(
    wm::AbstractDirectedModel, component_name::String; nw::Int=wm.cnw,
    bounded::Bool=true, report::Bool=true)
    # Store the corresponding component symbol.
    comp_sym = Symbol(component_name)

    # Initialize variables associated with positive flows.
    qp = var(wm, nw)[Symbol("qp_" * component_name)] = JuMP.@variable(
        wm.model, [a in ids(wm, nw, comp_sym)], lower_bound=0.0, base_name="$(nw)_qp",
        start=comp_start_value(ref(wm, nw, comp_sym, a), "qp_start"))

    # Initialize variables associated with negative flows.
    qn = var(wm, nw)[Symbol("qn_" * component_name)] = JuMP.@variable(
        wm.model, [a in ids(wm, nw, comp_sym)], lower_bound=0.0, base_name="$(nw)_qn",
        start=comp_start_value(ref(wm, nw, comp_sym, a), "qn_start"))

    if bounded # Bound flow-related variables if desired.
        q_lb, q_ub = calc_flow_bounds(wm, nw)

        for (a, comp) in ref(wm, nw, comp_sym)
            JuMP.set_upper_bound(qp[a], max(0.0, maximum(q_ub[component_name][a])))
            JuMP.set_upper_bound(qn[a], max(0.0, -minimum(q_lb[component_name][a])))
        end
    end

    # Report positive directed flow values as part of the solution.
    report && sol_component_value(wm, nw, comp_sym, :qp, ids(wm, nw, comp_sym), qp)

    # Report negative directed flow values as part of the solution.
    report && sol_component_value(wm, nw, comp_sym, :qn, ids(wm, nw, comp_sym), qn)

    # Create expressions capturing the relationships among q, qp, and qn.
    q = var(wm, nw)[Symbol("q_" * component_name)] = JuMP.@expression(
        wm.model, [a in ids(wm, nw, comp_sym)], qp[a] - qn[a])

    # Report flow expression values as part of the solution.
    report && sol_component_value(wm, nw, comp_sym, :q, ids(wm, nw, comp_sym), q)
end


"Create flow-related variables common to all directed flow models for edge-type components."
function variable_flow(wm::AbstractDirectedModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    for name in ["check_valve", "pipe", "pressure_reducing_valve", "pump", "shutoff_valve"]
        # Create directed flow (`qp` and `qn`) variables for each component.
        _variable_component_flow(wm, name; nw=nw, bounded=bounded, report=report)

        # Create directed head difference (`dhp` and `dhn`) variables for each component.
        _variable_component_head_difference(wm, name; nw=nw, bounded=bounded, report=report)

        # Create directed flow binary direction variables (`y`) for each component.
        _variable_component_direction(wm, name; nw=nw, report=report)
    end
end


"Create network design flow variables for directed flow formulations."
function variable_flow_des(wm::AbstractDirectedModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    # Create dictionary for undirected design flow variables (qp_des and qn_des).
    qp_des = var(wm, nw)[:qp_des] = Dict{Int,Array{JuMP.VariableRef}}()
    qn_des = var(wm, nw)[:qn_des] = Dict{Int,Array{JuMP.VariableRef}}()

    # Initialize the variables. (The default start value of _q_eps is crucial.)
    for a in ids(wm, nw, :pipe_des)
        var(wm, nw, :qp_des)[a] = JuMP.@variable(wm.model,
            [r in 1:length(ref(wm, nw, :resistance, a))], lower_bound=0.0,
            base_name="$(nw)_qp_des[$(a)]",
            start=comp_start_value(ref(wm, nw, :pipe_des, a), "qp_des_start", r, _q_eps))

        var(wm, nw, :qn_des)[a] = JuMP.@variable(wm.model,
            [r in 1:length(ref(wm, nw, :resistance, a))], lower_bound=0.0,
            base_name="$(nw)_qn_des[$(a)]",
            start=comp_start_value(ref(wm, nw, :pipe_des, a), "qn_des_start", r, _q_eps))
    end

    if bounded # If the variables are bounded, apply the bounds.
        q_lb, q_ub = calc_flow_bounds(wm, nw)

        for a in ids(wm, nw, :pipe_des)
            for r in 1:length(ref(wm, nw, :resistance, a))
                JuMP.set_upper_bound(qp_des[a][r], max(0.0, q_ub["pipe"][a][r]))
                JuMP.set_upper_bound(qn_des[a][r], max(0.0, -q_lb["pipe"][a][r]))
            end
        end
    end

    # Create expressions capturing the relationships among q, qp_des, and qn_des.
    q = var(wm, nw)[:q] = JuMP.@expression(
        wm.model, [a in ids(wm, nw, :pipe_des)],
        sum(var(wm, nw, :qp_des, a)) - sum(var(wm, nw, :qn_des, a)))

    # Initialize the solution reporting data structures.
    report && sol_component_value(wm, nw, :pipe_des, :q, ids(wm, nw, :pipe_des), q)

    # Create resistance binary variables.
    variable_resistance(wm, nw=nw)
end


function variable_pump_operation(wm::AbstractDirectedModel; nw::Int=wm.cnw, report::Bool=true)
    # If the number of breakpoints is not positive, return.
    pump_breakpoints = get(wm.ext, :pump_breakpoints, 0)

    if pump_breakpoints > 0
        # Create weights involved in convex combination constraints.
        lambda = var(wm, nw)[:lambda_pump] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pump), k in 1:pump_breakpoints],
            base_name="$(nw)_lambda", lower_bound=0.0, upper_bound=1.0,
            start=comp_start_value(ref(wm, nw, :pump, a), "lambda_start", k))

        # Create binary variables involved in convex combination constraints.
        x_pw = var(wm, nw)[:x_pw_pump] = JuMP.@variable(wm.model,
            [a in ids(wm, nw, :pump), k in 1:pump_breakpoints-1],
            base_name="$(nw)_x_pw", binary=true,
            start=comp_start_value(ref(wm, nw, :pump, a), "x_pw_start", k))
    end
end


function constraint_check_valve_common(wm::AbstractDirectedModel, n::Int, a::Int, node_fr::Int, node_to::Int)
    # Get common flow variables.
    qp, qn = var(wm, n, :qp_check_valve, a), var(wm, n, :qn_check_valve, a)
    y, z = var(wm, n, :y_check_valve, a), var(wm, n, :z_check_valve, a)

    # Fix flow in the negative direction to zero.
    JuMP.fix(qn, 0.0, force=true)

    # If the check valve is open, flow must be appreciably nonnegative.
    c_1 = JuMP.@constraint(wm.model, qp <= JuMP.upper_bound(qp) * z)
    c_2 = JuMP.@constraint(wm.model, qp >= _q_eps * z)

    # Get common head variables and associated data.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)
    dhp, dhn = var(wm, n, :dhp_check_valve, a), var(wm, n, :dhn_check_valve, a)
    dhp_ub, dhn_ub = JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)

    # When the check valve is closed, positive head loss is not possible.
    c_3 = JuMP.@constraint(wm.model, dhp <= dhp_ub * z)

    # When the check valve is open, negative head loss is not possible.
    c_4 = JuMP.@constraint(wm.model, dhn <= dhn_ub * (1.0 - z))

    # Constrain head differences based on direction.
    c_5 = JuMP.@constraint(wm.model, dhp <= dhp_ub * y)
    c_6 = JuMP.@constraint(wm.model, dhn <= dhn_ub * (1.0 - y))

    # Constrain directed flows based on direction.
    c_7 = JuMP.@constraint(wm.model, qp <= JuMP.upper_bound(qp) * y)

    # Relate head differences with head variables
    c_8 = JuMP.@constraint(wm.model, dhp - dhn == h_i - h_j)

    # Append the constraint array.
    append!(con(wm, n, :check_valve, a), [c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8])
end


function constraint_sv_common(wm::AbstractDirectedModel, n::Int, a::Int, node_fr::Int, node_to::Int)
    # Get flow and shutoff valve status variables.
    y, z = var(wm, n, :y_shutoff_valve, a), var(wm, n, :z_shutoff_valve, a)
    qp, qn = var(wm, n, :qp_shutoff_valve, a), var(wm, n, :qn_shutoff_valve, a)

    # Get common head variables and associated data.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)
    dhp, dhn = var(wm, n, :dhp_shutoff_valve, a), var(wm, n, :dhn_shutoff_valve, a)
    dhp_ub, dhn_ub = JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)

    # Constrain head differences based on direction.
    c_1 = JuMP.@constraint(wm.model, dhp <= dhp_ub * y)
    c_2 = JuMP.@constraint(wm.model, dhn <= dhn_ub * (1.0 - y))

    # Constrain directed flows based on direction.
    c_3 = JuMP.@constraint(wm.model, qp <= JuMP.upper_bound(qp) * y)
    c_4 = JuMP.@constraint(wm.model, qn <= JuMP.upper_bound(qn) * (1.0 - y))

    # Constrain directed flows based on shutoff valve status.
    c_5 = JuMP.@constraint(wm.model, qp <= JuMP.upper_bound(qp) * z)
    c_6 = JuMP.@constraint(wm.model, qn <= JuMP.upper_bound(qn) * z)

    # Relate head differences with head variables
    c_7 = JuMP.@constraint(wm.model, dhp - dhn == h_i - h_j)

    # Append the constraint array.
    append!(con(wm, n, :sv, a), [c_1, c_2, c_3, c_4, c_5, c_6, c_7])
end


function constraint_pipe_common(wm::AbstractDirectedModel, n::Int, a::Int, node_fr::Int, node_to::Int, alpha::Float64, L::Float64, r::Float64)
    # Get common flow variables and associated data.
    y = var(wm, n, :y_pipe, a)
    qp, qn = var(wm, n, :qp_pipe, a), var(wm, n, :qn_pipe, a)
    qp_ub, qn_ub = JuMP.upper_bound(qp), JuMP.upper_bound(qn)

    # Constrain directed flow variables based on direction.
    c_1 = JuMP.@constraint(wm.model, qp <= qp_ub * y)
    c_2 = JuMP.@constraint(wm.model, qn <= qn_ub * (1.0 - y))

    # Get common head variables and associated data.
    dhp, dhn = var(wm, n, :dhp_pipe, a), var(wm, n, :dhn_pipe, a)
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)
    dhp_ub, dhn_ub = JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)

    # Constrain directed head variables based on direction.
    c_3 = JuMP.@constraint(wm.model, dhp <= dhp_ub * y)
    c_4 = JuMP.@constraint(wm.model, dhn <= dhn_ub * (1.0 - y))
    c_5 = JuMP.@constraint(wm.model, dhp - dhn == h_i - h_j)

    # Append the constraint array.
    append!(con(wm, n, :pipe)[a], [c_1, c_2, c_3, c_4, c_5])
end


function constraint_prv_common(wm::AbstractDirectedModel, n::Int, a::Int, node_fr::Int, node_to::Int, h_prv::Float64)
    # Get common flow variables.
    qp = var(wm, n, :qp_pressure_reducing_valve, a)
    y = var(wm, n, :y_pressure_reducing_valve, a)
    z = var(wm, n, :z_pressure_reducing_valve, a)

    # If the pressure reducing valve is open, flow must be appreciably positive.
    c_1 = JuMP.@constraint(wm.model, qp <= JuMP.upper_bound(qp) * z)
    c_2 = JuMP.@constraint(wm.model, qp <= JuMP.upper_bound(qp) * y)
    c_3 = JuMP.@constraint(wm.model, qp >= _q_eps * z)

    # Get common head variables and associated data.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)
    dhp = var(wm, n, :dhp_pressure_reducing_valve, a)
    dhn = var(wm, n, :dhn_pressure_reducing_valve, a)
    dhp_ub, dhn_ub = JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)

    # When the pressure reducing valve is open, the head at node j is predefined.
    h_lb, h_ub = JuMP.lower_bound(h_j), JuMP.upper_bound(h_j)
    c_4 = JuMP.@constraint(wm.model, h_j >= (1.0 - z) * h_lb + z * h_prv)
    c_5 = JuMP.@constraint(wm.model, h_j <= (1.0 - z) * h_ub + z * h_prv)

    # Constrain directed head differences based on status.
    c_6 = JuMP.@constraint(wm.model, dhp <= dhp_ub * z)
    c_7 = JuMP.@constraint(wm.model, dhn <= dhn_ub * (1.0 - z))
    c_8 = JuMP.@constraint(wm.model, dhp - dhn == h_i - h_j)

    # Append the constraint array.
    append!(con(wm, n, :prv, a), [c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8])
end


function constraint_pump_common(wm::AbstractDirectedModel, n::Int, a::Int, node_fr::Int, node_to::Int, pc::Array{Float64})
    # Gather variables common to all formulations involving pumps.
    y, z = var(wm, n, :y_pump, a), var(wm, n, :z_pump, a)
    qp, g = var(wm, n, :qp_pump, a), var(wm, n, :g, a)
    dhp, dhn = var(wm, n, :dhp_pump, a), var(wm, n, :dhn_pump, a)

    # If the pump is off, flow across the pump must be zero.
    qp_ub = JuMP.upper_bound(qp)
    c_1 = JuMP.@constraint(wm.model, qp <= qp_ub * z)
    c_2 = JuMP.@constraint(wm.model, qp >= _q_eps * z)

    # If the pump is off, decouple the head difference relationship.
    dhp_ub, dhn_ub = JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)
    c_3 = JuMP.@constraint(wm.model, dhp <= dhp_ub * (1.0 - z))
    c_4 = JuMP.@constraint(wm.model, dhn <= g + dhn_ub * (1.0 - z))
    c_5 = JuMP.@constraint(wm.model, dhn >= g)

    # If the pump is on, flow across the pump must be nonnegative.
    c_6 = JuMP.@constraint(wm.model, y == z)

    # Equate head difference variables with heads.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)
    c_7 = JuMP.@constraint(wm.model, dhp - dhn == h_i - h_j)

    # Add a linear lower bound on the head gain approximation.
    g_1, g_2 = pc[3], pc[1]*qp_ub^2 + pc[2]*qp_ub + pc[3]
    g_lb_line = (g_2 - g_1) * inv(qp_ub) * qp + g_1 * z
    c_8 = JuMP.@constraint(wm.model, g_lb_line <= g)

    # Append the constraint array.
    append!(con(wm, n, :pump, a), [c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8])
end


function constraint_flow_direction_selection_des(wm::AbstractDirectedModel, n::Int, a::Int, pipe_resistances)
    y = var(wm, n, :y_pipe, a)

    for r_id in 1:length(pipe_resistances)
        qp, qn = var(wm, n, :qp_des, a)[r_id], var(wm, n, :qn_des, a)[r_id]
        c_p = JuMP.@constraint(wm.model, qp <= JuMP.upper_bound(qp) * y)
        c_n = JuMP.@constraint(wm.model, qn <= JuMP.upper_bound(qn) * (1.0 - y))
        append!(con(wm, n, :head_loss)[a], [c_p, c_n])
    end
end


function constraint_head_loss_ub_cv(wm::AbstractDirectedModel, n::Int, a::Int, alpha::Float64, L::Float64, r::Float64)
    qp, dhp = var(wm, n, :qp_check_valve, a), var(wm, n, :dhp_check_valve, a)
    rhs_p = r * JuMP.upper_bound(qp)^(alpha - 1.0) * qp
    c_p = JuMP.@constraint(wm.model, inv(L) * dhp <= rhs_p)
    append!(con(wm, n, :head_loss)[a], [c_p])
end


function constraint_shutoff_valve_head_loss_ub(wm::AbstractDirectedModel, n::Int, a::Int, alpha::Float64, L::Float64, r::Float64)
    z = var(wm, n, :z_shutoff_valve, a)

    qp, dhp = var(wm, n, :qp_shutoff_valve, a), var(wm, n, :dhp_shutoff_valve, a)
    rhs_p = JuMP.upper_bound(dhp) * (1.0 - z) + L*r*JuMP.upper_bound(qp)^(alpha - 1.0) * qp
    c_p = JuMP.@constraint(wm.model, dhp <= rhs_p)

    qn, dhn = var(wm, n, :qn_shutoff_valve, a), var(wm, n, :dhn_shutoff_valve, a)
    rhs_n = JuMP.upper_bound(dhn) * (1.0 - z) + L*r*JuMP.upper_bound(qn)^(alpha - 1.0) * qn
    c_n = JuMP.@constraint(wm.model, dhn <= rhs_n)

    append!(con(wm, n, :head_loss)[a], [c_p, c_n])
end


function constraint_pipe_head_loss_ub(wm::AbstractDirectedModel, n::Int, a::Int, alpha::Float64, L::Float64, r::Float64)
    qp, dhp = var(wm, n, :qp_pipe, a), var(wm, n, :dhp_pipe, a)
    rhs_p = r * JuMP.upper_bound(qp)^(alpha - 1.0) * qp
    c_p = JuMP.@constraint(wm.model, inv(L) * dhp <= rhs_p)

    qn, dhn = var(wm, n, :qn_pipe, a), var(wm, n, :dhn_pipe, a)
    rhs_n = r * JuMP.upper_bound(qn)^(alpha - 1.0) * qn
    c_n = JuMP.@constraint(wm.model, inv(L) * dhn <= rhs_n)

    append!(con(wm, n, :head_loss)[a], [c_p, c_n])
end


function constraint_pipe_head_loss_ub_des(wm::AbstractDirectedModel, n::Int, a::Int, alpha::Float64, L::Float64, pipe_resistances)
    dhp = var(wm, n, :dhp_pipe, a)
    qp_des = var(wm, n, :qp_des, a)
    qp_des_ub = JuMP.upper_bound.(qp_des)
    slopes_p = pipe_resistances .* qp_des_ub.^(alpha - 1.0)
    c_p = JuMP.@constraint(wm.model, inv(L)*dhp <= sum(slopes_p .* qp_des))

    dhn = var(wm, n, :dhn_pipe, a)
    qn_des = var(wm, n, :qn_des, a)
    qn_des_ub = JuMP.upper_bound.(qn_des)
    slopes_n = pipe_resistances .* qn_des_ub.^(alpha - 1.0)
    c_n = JuMP.@constraint(wm.model, inv(L)*dhn <= sum(slopes_n .* qn_des))

    append!(con(wm, n, :head_loss)[a], [c_p, c_n])
end


function constraint_resistance_selection_des(wm::AbstractDirectedModel, n::Int, a::Int, pipe_resistances)
    c = JuMP.@constraint(wm.model, sum(var(wm, n, :x_res, a)) == 1.0)
    append!(con(wm, n, :head_loss)[a], [c])

    for (r_id, r) in enumerate(pipe_resistances)
        x_res = var(wm, n, :x_res, a)[r_id]

        qp_des = var(wm, n, :qp_des, a)[r_id]
        qp_ub = JuMP.upper_bound(qp_des)
        c_p = JuMP.@constraint(wm.model, qp_des <= qp_ub * x_res)

        qn_des = var(wm, n, :qn_des, a)[r_id]
        qn_ub = JuMP.upper_bound(qn_des)
        c_n = JuMP.@constraint(wm.model, qn_des <= qn_ub * x_res)

        append!(con(wm, n, :head_loss)[a], [c_p, c_n])
    end
end


"Constraint to ensure at least one direction is set to take flow away from a source."
function constraint_source_directionality(
    wm::AbstractDirectedModel, n::Int, i::Int, check_valve_fr::Array{Int64,1},
    check_valve_to::Array{Int64,1}, pipe_fr::Array{Int64,1}, pipe_to::Array{Int64,1},
    pump_fr::Array{Int64,1}, pump_to::Array{Int64,1},
    pressure_reducing_valve_fr::Array{Int64,1}, pressure_reducing_valve_to::Array{Int64,1},
    shutoff_valve_fr::Array{Int64,1}, shutoff_valve_to::Array{Int64,1})
    # Collect direction variable references per component.
    y_check_valve = var(wm, n, :y_check_valve)
    y_pipe, y_pump = var(wm, n, :y_pipe), var(wm, n, :y_pump)
    y_pressure_reducing_valve = var(wm, n, :y_pressure_reducing_valve)
    y_shutoff_valve = var(wm, n, :y_shutoff_valve)

    y_out = JuMP.@expression(wm.model,
            sum(y_check_valve[a] for a in check_valve_fr) +
            sum(y_pipe[a] for a in pipe_fr) + sum(y_pump[a] for a in pump_fr) +
            sum(y_pressure_reducing_valve[a] for a in pressure_reducing_valve_fr) +
            sum(y_shutoff_valve[a] for a in shutoff_valve_fr))

    y_in = JuMP.@expression(wm.model,
            sum(y_check_valve[a] for a in check_valve_to) +
            sum(y_pipe[a] for a in pipe_to) + sum(y_pump[a] for a in pump_to) +
            sum(y_pressure_reducing_valve[a] for a in pressure_reducing_valve_to) +
            sum(y_shutoff_valve[a] for a in shutoff_valve_to))

    y_in_length = length(check_valve_to) + length(pipe_to) + length(pump_to) +
                  length(pressure_reducing_valve_to) + length(shutoff_valve_to)

    # Add the source flow direction constraint.
    c = JuMP.@constraint(wm.model, y_out - y_in >= 1.0 - y_in_length)
    con(wm, n, :source_flow)[i] = c
end


"Constraint to ensure at least one direction is set to take flow to a junction with demand."
function constraint_sink_directionality(
    wm::AbstractDirectedModel, n::Int, i::Int, check_valve_fr::Array{Int64,1},
    check_valve_to::Array{Int64,1}, pipe_fr::Array{Int64,1}, pipe_to::Array{Int64,1},
    pump_fr::Array{Int64,1}, pump_to::Array{Int64,1},
    pressure_reducing_valve_fr::Array{Int64,1}, pressure_reducing_valve_to::Array{Int64,1},
    shutoff_valve_fr::Array{Int64,1}, shutoff_valve_to::Array{Int64,1})
    # Collect direction variable references per component.
    y_check_valve = var(wm, n, :y_check_valve)
    y_pipe, y_pump = var(wm, n, :y_pipe), var(wm, n, :y_pump)
    y_pressure_reducing_valve = var(wm, n, :y_pressure_reducing_valve)
    y_shutoff_valve = var(wm, n, :y_shutoff_valve)

    y_out = JuMP.@expression(wm.model,
            sum(y_check_valve[a] for a in check_valve_fr) +
            sum(y_pipe[a] for a in pipe_fr) + sum(y_pump[a] for a in pump_fr) +
            sum(y_pressure_reducing_valve[a] for a in pressure_reducing_valve_fr) +
            sum(y_shutoff_valve[a] for a in shutoff_valve_fr))

    y_in = JuMP.@expression(wm.model,
          sum(y_check_valve[a] for a in check_valve_to) +
          sum(y_pipe[a] for a in pipe_to) + sum(y_pump[a] for a in pump_to) +
          sum(y_pressure_reducing_valve[a] for a in pressure_reducing_valve_to) +
          sum(y_shutoff_valve[a] for a in shutoff_valve_to))

    y_out_length = length(check_valve_fr) + length(pipe_fr) + length(pump_fr) +
                   length(pressure_reducing_valve_fr) + length(shutoff_valve_fr)

    # Add the sink flow direction constraint.
    c = JuMP.@constraint(wm.model, y_in - y_out >= 1.0 - y_out_length)
    con(wm, n, :sink_directionality)[i] = c
end
