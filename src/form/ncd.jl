# Define common NCD (nonlinear nonconvex and direction-based) implementations
# of water distribution network constraints, which use directed flow variables.

# Constraints and variables common to all formulations with directed flows. In these
# formulations, the variables qp correspond to flow from i to j, and the variables qn
# correspond to flow from j to i. That is, when qp is nonzero, qn should be zero, and when
# qn is nonzero, qp should be zero.


"Initialize variables associated with flow direction. If this variable is equal to one, the
flow direction is from i to j. If it is equal to zero, the flow direction is from j to i."
function _variable_component_direction(
    wm::AbstractNCDModel,
    component_name::String;
    nw::Int = nw_id_default,
    report::Bool = true,
)
    # Store the corresponding component symbol.
    comp_sym = Symbol(component_name)

    # Initialize variables associated with positive flows.
    y =
        var(wm, nw)[Symbol("y_" * component_name)] = JuMP.@variable(
            wm.model,
            [a in ids(wm, nw, comp_sym)],
            binary = true,
            base_name = "$(nw)_y",
            start = comp_start_value(
                ref(wm, nw, comp_sym, a),
                "y_$(component_name)_start",
                1.0,
            )
        )

    for (a, comp) in ref(wm, nw, comp_sym)
        _fix_indicator_variable(y[a], comp, "y")
    end

    # Report back flow values as part of the solution.
    report && sol_component_value(wm, nw, comp_sym, :y, ids(wm, nw, comp_sym), y)
end


"Create head differences variables common to all directed flow models for a component."
function _variable_component_head_difference(
    wm::AbstractNCDModel,
    component_name::String;
    nw::Int = nw_id_default,
    bounded::Bool = true,
    report::Bool = true,
)
    # Store the corresponding component symbol.
    comp_sym = Symbol(component_name)

    # Initialize directed variables associated with positive flows.
    dhp =
        var(wm, nw)[Symbol("dhp_" * component_name)] = JuMP.@variable(
            wm.model,
            [a in ids(wm, nw, comp_sym)],
            lower_bound = 0.0,
            base_name = "$(nw)_dhp",
            start = comp_start_value(ref(wm, nw, comp_sym, a), "dhp_start")
        )

    # Initialize directed variables associated with negative flows.
    dhn =
        var(wm, nw)[Symbol("dhn_" * component_name)] = JuMP.@variable(
            wm.model,
            [a in ids(wm, nw, comp_sym)],
            lower_bound = 0.0,
            base_name = "$(nw)_dhn",
            start = comp_start_value(ref(wm, nw, comp_sym, a), "dhn_start")
        )

    if bounded # Bound flow-related variables if desired.
        # Set lower and upper bounds on head differences.
        for (a, comp) in ref(wm, nw, comp_sym)
            node_fr = ref(wm, nw, :node, comp["node_fr"])
            node_to = ref(wm, nw, :node, comp["node_to"])

            dhp_max = max(0.0, node_fr["head_max"] - node_to["head_min"])
            JuMP.set_upper_bound(dhp[a], dhp_max)
            dhp_start = comp_start_value(comp, "dhp_start", 0.5 * dhp_max)
            JuMP.set_start_value(dhp[a], dhp_start)

            dhn_max = max(0.0, node_to["head_max"] - node_fr["head_min"])
            JuMP.set_upper_bound(dhn[a], dhn_max)
            dhn_start = comp_start_value(comp, "dhn_start", 0.0)
            JuMP.set_start_value(dhn[a], dhn_start)
        end
    end

    # Report positive head difference values as part of the solution.
    report && sol_component_value(wm, nw, comp_sym, :dhp, ids(wm, nw, comp_sym), dhp)

    # Report positive head difference values as part of the solution.
    report && sol_component_value(wm, nw, comp_sym, :dhn, ids(wm, nw, comp_sym), dhn)
end


"Create flow variables that are common to all directed flow models for a component."
function _variable_component_flow(
    wm::AbstractNCDModel,
    component_name::String;
    nw::Int = nw_id_default,
    bounded::Bool = true,
    report::Bool = true,
)
    # Store the corresponding component symbol.
    comp_sym = Symbol(component_name)

    # Get the transformed minimum flow.
    wm_data = get_wm_data(wm.data)
    flow_transform = _calc_flow_per_unit_transform(wm_data)
    flow_min_scaled = flow_transform(_FLOW_MIN)

    # Initialize variables associated with positive flows.
    qp =
        var(wm, nw)[Symbol("qp_" * component_name)] = JuMP.@variable(
            wm.model,
            [a in ids(wm, nw, comp_sym)],
            lower_bound = 0.0,
            base_name = "$(nw)_qp",
            start =
                comp_start_value(ref(wm, nw, comp_sym, a), "qp_start", flow_min_scaled)
        )

    # Initialize variables associated with negative flows.
    qn =
        var(wm, nw)[Symbol("qn_" * component_name)] = JuMP.@variable(
            wm.model,
            [a in ids(wm, nw, comp_sym)],
            lower_bound = 0.0,
            base_name = "$(nw)_qn",
            start =
                comp_start_value(ref(wm, nw, comp_sym, a), "qn_start", flow_min_scaled)
        )

    if bounded # Bound flow-related variables if desired.
        for (a, comp) in ref(wm, nw, comp_sym)
            qp_min = max(0.0, comp["flow_min"])
            JuMP.set_lower_bound(qp[a], qp_min)
            qp_max = max(0.0, comp["flow_max"])
            JuMP.set_upper_bound(qp[a], qp_max)
            qp_start = comp_start_value(comp, "qp_start", 0.5 * qp_max)
            JuMP.set_start_value(qp[a], qp_start)

            qn_min = max(0.0, -comp["flow_max"])
            JuMP.set_lower_bound(qn[a], qn_min)
            qn_max = max(0.0, -comp["flow_min"])
            JuMP.set_upper_bound(qn[a], qn_max)
            qn_start = comp_start_value(comp, "qn_start", 0.0)
            JuMP.set_start_value(qn[a], qn_start)
        end
    end

    # Report positive directed flow values as part of the solution.
    report && sol_component_value(wm, nw, comp_sym, :qp, ids(wm, nw, comp_sym), qp)

    # Report negative directed flow values as part of the solution.
    report && sol_component_value(wm, nw, comp_sym, :qn, ids(wm, nw, comp_sym), qn)

    # Create expressions capturing the relationships among q, qp, and qn.
    q =
        var(wm, nw)[Symbol("q_" * component_name)] =
            JuMP.@expression(wm.model, [a in ids(wm, nw, comp_sym)], qp[a] - qn[a])

    # Report flow expression values as part of the solution.
    report && sol_component_value(wm, nw, comp_sym, :q, ids(wm, nw, comp_sym), q)
end


"Create flow-related variables common to all directed flow models for node-connecting components."
function variable_flow(
    wm::AbstractNCDModel;
    nw::Int = nw_id_default,
    bounded::Bool = true,
    report::Bool = true,
)
    for name in _LINK_COMPONENTS
        # Create directed flow (`qp` and `qn`) variables for each component.
        _variable_component_flow(wm, name; nw = nw, bounded = bounded, report = report)

        # Create directed flow binary direction variables (`y`) for each component.
        _variable_component_direction(wm, name; nw = nw, report = report)
    end

    for name in ["des_pipe", "pipe"]
        # Create directed head difference (`dhp` and `dhn`) variables for each component.
        _variable_component_head_difference(
            wm,
            name;
            nw = nw,
            bounded = bounded,
            report = report,
        )
    end
end


function constraint_pipe_flow(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # Get common flow variables and associated data.
    qp, qn, y = var(wm, n, :qp_pipe, a), var(wm, n, :qn_pipe, a), var(wm, n, :y_pipe, a)

    # Constrain directed flow variables based on direction.
    qp_ub, qn_ub = JuMP.upper_bound(qp), JuMP.upper_bound(qn)
    c_1 = JuMP.@constraint(wm.model, qp <= qp_ub * y)
    c_2 = JuMP.@constraint(wm.model, qn <= qn_ub * (1.0 - y))

    # Add additional constraints based on active flows.
    qp_min_forward, qn_min_forward = max(0.0, q_min_forward), max(0.0, -q_max_reverse)
    c_3 = JuMP.@constraint(wm.model, qp >= qp_min_forward * y)
    c_4 = JuMP.@constraint(wm.model, qn >= qn_min_forward * (1.0 - y))

    # Append the :pipe_flow constraint array.
    append!(con(wm, n, :pipe_flow, a), [c_1, c_2, c_3, c_4])
end


function constraint_pipe_head(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
)
    # Get common flow variables and associated data.
    dhp, dhn, y = var(wm, n, :dhp_pipe, a), var(wm, n, :dhn_pipe, a), var(wm, n, :y_pipe, a)

    # Get head variables for from and to nodes.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # For pipes, the differences must satisfy lower and upper bounds.
    dhp_ub, dhn_ub = JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)
    c_1 = JuMP.@constraint(wm.model, dhp <= dhp_ub * y)
    c_2 = JuMP.@constraint(wm.model, dhn <= dhn_ub * (1.0 - y))
    c_3 = JuMP.@constraint(wm.model, dhp - dhn == h_i - h_j)

    # Append the :pipe_head constraint array.
    append!(con(wm, n, :pipe_head, a), [c_1, c_2, c_3])
end


function constraint_pipe_head_loss(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    exponent::Float64,
    L::Float64,
    r::Float64,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # Add constraints for positive flow and head difference.
    qp, dhp = var(wm, n, :qp_pipe, a), var(wm, n, :dhp_pipe, a)
    c_1 = JuMP.@NLconstraint(wm.model, r * head_loss(qp) <= dhp / L)
    c_2 = JuMP.@NLconstraint(wm.model, r * head_loss(qp) >= dhp / L)

    # Add constraints for negative flow and head difference.
    qn, dhn = var(wm, n, :qn_pipe, a), var(wm, n, :dhn_pipe, a)
    c_3 = JuMP.@NLconstraint(wm.model, r * head_loss(qn) <= dhn / L)
    c_4 = JuMP.@NLconstraint(wm.model, r * head_loss(qn) >= dhn / L)

    # Append the :pipe_head_loss constraint array.
    append!(con(wm, n, :pipe_head_loss)[a], [c_1, c_2, c_3, c_4])
end


function constraint_on_off_des_pipe_flow(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # Get des_pipe status variable.
    qp, qn = var(wm, n, :qp_des_pipe, a), var(wm, n, :qn_des_pipe, a)
    y, z = var(wm, n, :y_des_pipe, a), var(wm, n, :z_des_pipe, a)

    # If the des_pipe is inactive, flow must be zero.
    qp_ub, qn_ub = JuMP.upper_bound(qp), JuMP.upper_bound(qn)
    c_1 = JuMP.@constraint(wm.model, qp <= qp_ub * y)
    c_2 = JuMP.@constraint(wm.model, qn <= qn_ub * (1.0 - y))
    c_3 = JuMP.@constraint(wm.model, qp <= qp_ub * z)
    c_4 = JuMP.@constraint(wm.model, qn <= qn_ub * z)

    # Append the constraint array.
    append!(con(wm, n, :on_off_des_pipe_flow, a), [c_1, c_2, c_3, c_4])
end


function constraint_on_off_des_pipe_head(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
)
    # Get head difference variables for the des_pipe.
    dh = var(wm, n, :dh_des_pipe, a)
    dhp = var(wm, n, :dhp_des_pipe, a)
    dhn = var(wm, n, :dhn_des_pipe, a)

    # Get des_pipe direction and status variable.
    y = var(wm, n, :y_des_pipe, a)
    z = var(wm, n, :z_des_pipe, a)

    # If the des_pipe is off, decouple the head difference relationship.
    dhp_ub, dhn_ub = JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)
    c_1 = JuMP.@constraint(wm.model, dhp <= dhp_ub * y)
    c_2 = JuMP.@constraint(wm.model, dhn <= dhn_ub * (1.0 - y))
    c_3 = JuMP.@constraint(wm.model, dhp <= dhp_ub * z)
    c_4 = JuMP.@constraint(wm.model, dhn <= dhn_ub * z)
    c_5 = JuMP.@constraint(wm.model, dh == dhp - dhn)

    # Append the constraint array.
    append!(con(wm, n, :on_off_des_pipe_head, a), [c_1, c_2, c_3, c_4, c_5])
end


"Adds head loss constraint for a design pipe in the `NC` formulation."
function constraint_on_off_des_pipe_head_loss(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    exponent::Float64,
    L::Float64,
    r::Float64,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # Get flow and design status variables.
    qp, qn = var(wm, n, :qp_des_pipe, a), var(wm, n, :qn_des_pipe, a)
    dhp, dhn = var(wm, n, :dhp_des_pipe, a), var(wm, n, :dhn_des_pipe, a)

    # Add nonconvex constraint for the head loss relationship.
    c_1 = JuMP.@NLconstraint(wm.model, r * head_loss(qp) <= dhp / L)
    c_2 = JuMP.@NLconstraint(wm.model, r * head_loss(qp) >= dhp / L)
    c_3 = JuMP.@NLconstraint(wm.model, r * head_loss(qn) <= dhn / L)
    c_4 = JuMP.@NLconstraint(wm.model, r * head_loss(qn) >= dhn / L)

    # Append the :pipe_head_loss constraint array.
    append!(con(wm, n, :on_off_des_pipe_head_loss)[a], [c_1, c_2, c_3, c_4])
end


function constraint_des_pipe_flow(
    wm::AbstractNCDModel,
    n::Int,
    k::Int,
    node_fr::Int,
    node_to::Int,
    des_pipes::Vector{Int},
)
    y_des_pipe = var(wm, n, :y_des_pipe)
    lhs = sum(y_des_pipe[a] for a in des_pipes)
    rhs = length(des_pipes) * y_des_pipe[des_pipes[1]]
    c = JuMP.@constraint(wm.model, lhs == rhs) # All directions are the same.
    append!(con(wm, n, :des_pipe_flow)[k], [c])
end


function constraint_des_pipe_head(
    wm::AbstractNCDModel,
    n::Int,
    k::Int,
    node_fr::Int,
    node_to::Int,
    des_pipes::Vector{Int}
)
    # Collect relevant design pipe variables and summations.
    h_i, h_j = var(wm, n, :h, node_fr), var(wm, n, :h, node_to)

    # Relate directed head differences with actual head difference.
    dhp_sum = sum(var(wm, n, :dhp_des_pipe, a) for a in des_pipes)
    dhn_sum = sum(var(wm, n, :dhn_des_pipe, a) for a in des_pipes)
    c = JuMP.@constraint(wm.model, dhp_sum - dhn_sum == h_i - h_j)
    append!(con(wm, n, :des_pipe_head)[k], [c])
end


function constraint_on_off_pump_flow(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    q_min_forward::Float64,
)
    if(n==1)
        println("Running ncd pump flow constraint")
    end
    # Get pump status variable.
    qp = var(wm, n, :qp_pump, a)
    y = var(wm, n, :y_pump, a)
    z = var(wm, n, :z_pump, a)

    # If the pump is inactive, flow must be zero.
    qp_lb, qp_ub = q_min_forward, JuMP.upper_bound(qp)
    c_1 = JuMP.@constraint(wm.model, qp >= qp_lb * z)
    c_2 = JuMP.@constraint(wm.model, qp <= qp_ub * z)

    # If the pump is on, the flow direction must be positive.
    c_3 = JuMP.@constraint(wm.model, y >= z)

    # Append the constraint array.
    append!(con(wm, n, :on_off_pump_flow, a), [c_1, c_2, c_3])
end

function constraint_on_off_pump_flow_ne(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    q_min_forward::Float64,
)
    if(n == 1)
        println("Running NCD pump flow ne constraint")
    end
    # Get pump status variable.
    qp = var(wm, n, :qp_ne_pump, a)
    y = var(wm, n, :y_ne_pump, a)
    z = var(wm, n, :z_ne_pump, a)

    # If the pump is inactive, flow must be zero.
    qp_lb, qp_ub = q_min_forward, JuMP.upper_bound(qp)
    c_1 = JuMP.@constraint(wm.model, qp >= qp_lb * z)
    c_2 = JuMP.@constraint(wm.model, qp <= qp_ub * z)

    # If the pump is on, the flow direction must be positive.
    c_3 = JuMP.@constraint(wm.model, y >= z)

    # Append the constraint array.
    append!(con(wm, n, :on_off_pump_flow_ne, a), [c_1, c_2, c_3])
end


function constraint_on_off_pump_head_gain(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    q_min_forward::Float64,
)
    # Gather pump flow, head gain, and status variables.
    qp = var(wm, n, :qp_pump, a)
    g = var(wm, n, :g_pump, a)
    z = var(wm, n, :z_pump, a)

    # Define the (relaxed) head gain relationship for the pump.
    head_curve_func = _calc_head_curve_function(ref(wm, n, :pump, a), z)
    c_1 = JuMP.@constraint(wm.model, head_curve_func(qp) <= g)
    c_2 = JuMP.@constraint(wm.model, head_curve_func(qp) >= g)
    append!(con(wm, n, :on_off_pump_head_gain, a), [c_1, c_2])
end

function constraint_on_off_pump_head_gain_ne(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    q_min_forward::Float64,
)
    if(n==1)
        println("Running NCD pump head gain ne constraint")
    end
    # Gather pump flow, head gain, and status variables.
    qp = var(wm, n, :qp_ne_pump, a)
    g = var(wm, n, :g_ne_pump, a)
    z = var(wm, n, :z_ne_pump, a)

    # Define the (relaxed) head gain relationship for the pump.
    head_curve_func = _calc_head_curve_function(ref(wm, n, :ne_pump, a), z)
    c_1 = JuMP.@constraint(wm.model, head_curve_func(qp) <= g)
    c_2 = JuMP.@constraint(wm.model, head_curve_func(qp) >= g)
    append!(con(wm, n, :on_off_pump_head_gain_ne, a), [c_1, c_2])
end

function constraint_on_off_pump_power(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    q_min_forward::Float64,
)
    # if(n==1)
        println("Running ncd pump $a power constraint")
    # end
    # Gather pump flow, power, and status variables.
    q = var(wm, n, :qp_pump, a)
    P = var(wm, n, :P_pump, a)
    z = var(wm, n, :z_pump, a)

    # Add constraint equating power with respect to the power curve.
    power_qa = _calc_pump_power_quadratic_approximation(wm, n, a, z)
    c_1 = JuMP.@constraint(wm.model, power_qa(q) <= P)
    c_2 = JuMP.@constraint(wm.model, power_qa(q) >= P)

    # Append the :on_off_pump_power constraint array.
    append!(con(wm, n, :on_off_pump_power)[a], [c_1, c_2])
end

function constraint_on_off_pump_power_ne(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    q_min_forward::Float64,
)
    if(n == 1)
        println("Running NCD pump power ne constraint")
    end
    # Gather pump flow, power, and status variables.
    q = var(wm, n, :qp_ne_pump, a)
    P = var(wm, n, :P_ne_pump, a)
    z = var(wm, n, :z_ne_pump, a)

    # Add constraint equating power with respect to the power curve.
    power_qa = _calc_pump_power_quadratic_approximation_ne(wm, n, a, z)
    c_1 = JuMP.@constraint(wm.model, power_qa(q) <= P)
    c_2 = JuMP.@constraint(wm.model, power_qa(q) >= P)

    # Append the :on_off_pump_power constraint array.
    append!(con(wm, n, :on_off_pump_power_ne)[a], [c_1, c_2])
end


function constraint_on_off_regulator_flow(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    q_min_forward::Float64,
)
    # Get regulator flow, status, and direction variables.
    qp, z = var(wm, n, :qp_regulator, a), var(wm, n, :z_regulator, a)
    y = var(wm, n, :y_regulator, a) # Regulator flow direction.

    # If the regulator is closed, flow must be zero.
    qp_lb, qp_ub = max(JuMP.lower_bound(qp), q_min_forward), JuMP.upper_bound(qp)
    c_1 = JuMP.@constraint(wm.model, qp >= qp_lb * z)
    c_2 = JuMP.@constraint(wm.model, qp <= qp_ub * z)
    c_3 = JuMP.@constraint(wm.model, qp <= qp_ub * y)

    # Append the :on_off_regulator_flow constraint array.
    append!(con(wm, n, :on_off_regulator_flow, a), [c_1, c_2, c_3])
end


function constraint_on_off_valve_flow(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # Get valve flow, direction, and status variables.
    qp, qn = var(wm, n, :qp_valve, a), var(wm, n, :qn_valve, a)
    y, z = var(wm, n, :y_valve, a), var(wm, n, :z_valve, a)

    # The valve flow is constrained by direction and status.
    qp_ub, qn_ub = JuMP.upper_bound(qp), JuMP.upper_bound(qn)
    c_1 = JuMP.@constraint(wm.model, qp <= qp_ub * y)
    c_2 = JuMP.@constraint(wm.model, qn <= qn_ub * (1.0 - y))
    c_3 = JuMP.@constraint(wm.model, qp <= qp_ub * z)
    c_4 = JuMP.@constraint(wm.model, qn <= qn_ub * z)

    # Add additional constraints based on active flows.
    qp_min_forward, qn_min_forward = max(0.0, q_min_forward), max(0.0, -q_max_reverse)
    c_5 = JuMP.@constraint(wm.model, qp >= qp_min_forward * (y + z - 1.0))
    c_6 = JuMP.@constraint(wm.model, qn >= qn_min_forward * (z - y))

    # Append the constraint array.
    append!(con(wm, n, :on_off_valve_flow, a), [c_1, c_2, c_3, c_4, c_5, c_6])
end


function constraint_short_pipe_flow(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # Get short pipe flow and direction variables.
    qp, qn = var(wm, n, :qp_short_pipe, a), var(wm, n, :qn_short_pipe, a)
    y = var(wm, n, :y_short_pipe, a) # Binary direction variable.

    # The valve flow is constrained by direction and status.
    qp_ub, qn_ub = JuMP.upper_bound(qp), JuMP.upper_bound(qn)
    c_1 = JuMP.@constraint(wm.model, qp <= qp_ub * y)
    c_2 = JuMP.@constraint(wm.model, qn <= qn_ub * (1.0 - y))

    # Add additional constraints based on active flows.
    qp_min_forward, qn_min_forward = max(0.0, q_min_forward), max(0.0, -q_max_reverse)
    c_3 = JuMP.@constraint(wm.model, qp >= qp_min_forward * y)
    c_4 = JuMP.@constraint(wm.model, qn >= qn_min_forward * (1.0 - y))

    # Append the :short_pipe_flow constraint array.
    append!(con(wm, n, :short_pipe_flow, a), [c_1, c_2, c_3, c_4])
end


function constraint_short_pipe_flow_ne(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
    # Get expansion short pipe flow, direction, and status variables.
    qp, qn = var(wm, n, :qp_ne_short_pipe, a), var(wm, n, :qn_ne_short_pipe, a)
    y, z = var(wm, n, :y_ne_short_pipe, a), var(wm, n, :z_ne_short_pipe, a)

    # The expansion short pipe's flow is constrained by direction and build status.
    qp_ub, qn_ub = JuMP.upper_bound(qp), JuMP.upper_bound(qn)
    c_1 = JuMP.@constraint(wm.model, qp <= qp_ub * y)
    c_2 = JuMP.@constraint(wm.model, qn <= qn_ub * (1.0 - y))
    c_3 = JuMP.@constraint(wm.model, qp <= qp_ub * z)
    c_4 = JuMP.@constraint(wm.model, qn <= qn_ub * z)

    # Add additional constraints based on active flows.
    qp_min_forward, qn_min_forward = max(0.0, q_min_forward), max(0.0, -q_max_reverse)
    c_5 = JuMP.@constraint(wm.model, qp >= qp_min_forward * (y + z - 1.0))
    c_6 = JuMP.@constraint(wm.model, qn >= qn_min_forward * (z - y))

    # Append the constraint array.
    append!(con(wm, n, :short_pipe_flow_ne, a), [c_1, c_2, c_3, c_4, c_5, c_6])
end


function _gather_directionality_data(
    wm::AbstractNCDModel,
    n::Int,
    pipe_fr::Vector{Int},
    pipe_to::Vector{Int},
    des_pipe_fr::Vector{Int},
    des_pipe_to::Vector{Int},
    pump_fr::Vector{Int},
    pump_to::Vector{Int},
    ne_pump_fr::Vector{Int},
    ne_pump_to::Vector{Int},
    regulator_fr::Vector{Int},
    regulator_to::Vector{Int},
    short_pipe_fr::Vector{Int},
    short_pipe_to::Vector{Int},
    ne_short_pipe_fr::Vector{Int},
    ne_short_pipe_to::Vector{Int},
    valve_fr::Vector{Int},
    valve_to::Vector{Int},
)
    # Collect direction variable references per component.
    y_pipe, y_des_pipe = var(wm, n, :y_pipe), var(wm, n, :y_des_pipe)
    y_pump, y_regulator = var(wm, n, :y_pump), var(wm, n, :y_regulator)
    y_ne_pump = var(wm, n, :y_ne_pump)
    y_short_pipe, y_ne_short_pipe = var(wm, n, :y_short_pipe), var(wm, n, :y_ne_short_pipe)
    y_valve = var(wm, n, :y_valve)

    sum_in = JuMP.@expression(
        wm.model,
        sum(y_pipe[a] for a in pipe_to) +
        sum(y_des_pipe[a] for a in des_pipe_to) +
        sum(y_pump[a] for a in pump_to) +
        sum(y_ne_pump[a] for a in ne_pump_to) +
        sum(y_regulator[a] for a in regulator_to) +
        sum(y_short_pipe[a] for a in short_pipe_to) +
        sum(y_ne_short_pipe[a] for a in ne_short_pipe_to) +
        sum(y_valve[a] for a in valve_to)
    )

    sum_out = JuMP.@expression(
        wm.model,
        sum(y_pipe[a] for a in pipe_fr) +
        sum(y_des_pipe[a] for a in des_pipe_fr) +
        sum(y_pump[a] for a in pump_fr) +
        sum(y_ne_pump[a] for a in ne_pump_fr) +
        sum(y_regulator[a] for a in regulator_fr) +
        sum(y_short_pipe[a] for a in short_pipe_fr) +
        sum(y_ne_short_pipe[a] for a in ne_short_pipe_fr) +
        sum(y_valve[a] for a in valve_fr)
    )

    # Get the in degree of node `i`.
    in_length =
        length(pipe_to) +
        length(des_pipe_to) +
        length(pump_to) +
        length(ne_pump_to) +
        length(regulator_to) +
        length(short_pipe_to) +
        length(ne_short_pipe_to) +
        length(valve_to)

    # Get the out degree of node `i`.
    out_length =
        length(pipe_fr) +
        length(des_pipe_fr) +
        length(pump_fr) +
        length(ne_pump_fr) +
        length(regulator_fr) +
        length(short_pipe_fr) +
        length(ne_short_pipe_fr) +
        length(valve_fr)

    return sum_in, sum_out, in_length, out_length
end


"Constraint to ensure at least one direction is set to take flow away from a source."
function constraint_intermediate_directionality(
    wm::AbstractNCDModel,
    n::Int,
    i::Int,
    pipe_fr::Vector{Int},
    pipe_to::Vector{Int},
    des_pipe_fr::Vector{Int},
    des_pipe_to::Vector{Int},
    pump_fr::Vector{Int},
    pump_to::Vector{Int},
    ne_pump_fr::Vector{Int},
    ne_pump_to::Vector{Int},
    regulator_fr::Vector{Int},
    regulator_to::Vector{Int},
    short_pipe_fr::Vector{Int},
    short_pipe_to::Vector{Int},
    ne_short_pipe_fr::Vector{Int},
    ne_short_pipe_to::Vector{Int},
    valve_fr::Vector{Int},
    valve_to::Vector{Int},
)
    # Gather data required to build the constraint.
    sum_in, sum_out, in_length, out_length = _gather_directionality_data(
        wm,
        n,
        pipe_fr,
        pipe_to,
        des_pipe_fr,
        des_pipe_to,
        pump_fr,
        pump_to,
        ne_pump_fr,
        ne_pump_to,
        regulator_fr,
        regulator_to,
        short_pipe_fr,
        short_pipe_to,
        ne_short_pipe_fr,
        ne_short_pipe_to,
        valve_fr,
        valve_to,
    )

    # Add the node directionality constraint.
    if out_length == 1 && in_length == 1
        c = JuMP.@constraint(wm.model, sum_out - sum_in == 0.0)
        con(wm, n, :node_directionality)[i] = c
    elseif in_length + out_length == 2 && in_length * out_length == 0
        c = JuMP.@constraint(wm.model, sum_out + sum_in == 1.0)
        con(wm, n, :node_directionality)[i] = c
    end
end


"Constraint to ensure at least one direction is set to take flow away from a source."
function constraint_source_directionality(
    wm::AbstractNCDModel,
    n::Int,
    i::Int,
    pipe_fr::Vector{Int},
    pipe_to::Vector{Int},
    des_pipe_fr::Vector{Int},
    des_pipe_to::Vector{Int},
    pump_fr::Vector{Int},
    pump_to::Vector{Int},
    ne_pump_fr::Vector{Int},
    ne_pump_to::Vector{Int},
    regulator_fr::Vector{Int},
    regulator_to::Vector{Int},
    short_pipe_fr::Vector{Int},
    short_pipe_to::Vector{Int},
    ne_short_pipe_fr::Vector{Int},
    ne_short_pipe_to::Vector{Int},
    valve_fr::Vector{Int},
    valve_to::Vector{Int},
)
    # Gather data required to build the constraint.
    sum_in, sum_out, in_length, out_length = _gather_directionality_data(
        wm,
        n,
        pipe_fr,
        pipe_to,
        des_pipe_fr,
        des_pipe_to,
        pump_fr,
        pump_to,
        ne_pump_fr,
        ne_pump_to,
        regulator_fr,
        regulator_to,
        short_pipe_fr,
        short_pipe_to,
        ne_short_pipe_fr,
        ne_short_pipe_to,
        valve_fr,
        valve_to,
    )

    # Add the source flow direction constraint.
    c = JuMP.@constraint(wm.model, sum_out - sum_in >= 1.0 - in_length)
    con(wm, n, :node_directionality)[i] = c
end


"Constraint to ensure at least one direction is set to take flow to a node with demand."
function constraint_sink_directionality(
    wm::AbstractNCDModel,
    n::Int,
    i::Int,
    pipe_fr::Vector{Int},
    pipe_to::Vector{Int},
    des_pipe_fr::Vector{Int},
    des_pipe_to::Vector{Int},
    pump_fr::Vector{Int},
    pump_to::Vector{Int},
    ne_pump_fr::Vector{Int},
    ne_pump_to::Vector{Int},
    regulator_fr::Vector{Int},
    regulator_to::Vector{Int},
    short_pipe_fr::Vector{Int},
    short_pipe_to::Vector{Int},
    ne_short_pipe_fr::Vector{Int},
    ne_short_pipe_to::Vector{Int},
    valve_fr::Vector{Int},
    valve_to::Vector{Int},
)
    # Gather data required to build the constraint.
    sum_in, sum_out, in_length, out_length = _gather_directionality_data(
        wm,
        n,
        pipe_fr,
        pipe_to,
        des_pipe_fr,
        des_pipe_to,
        pump_fr,
        pump_to,
        ne_pump_fr,
        ne_pump_to,
        regulator_fr,
        regulator_to,
        short_pipe_fr,
        short_pipe_to,
        ne_short_pipe_fr,
        ne_short_pipe_to,
        valve_fr,
        valve_to,
    )

    # Add the sink flow direction constraint.
    c = JuMP.@constraint(wm.model, sum_in - sum_out >= 1.0 - out_length)
    con(wm, n, :node_directionality)[i] = c
end
