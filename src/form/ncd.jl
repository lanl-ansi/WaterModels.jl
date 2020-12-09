# Define common NCD (nonlinear nonconvex and direction-based) implementations
# of water distribution network constraints, which use directed flow variables.


function constraint_pipe_head_loss(
    wm::AbstractNCDModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Add constraints for positive flow and head difference.
    qp, dhp = var(wm, n, :qp_pipe, a), var(wm, n, :dhp_pipe, a)
    c_1 = JuMP.@NLconstraint(wm.model, r * head_loss(qp) <= inv(L) * dhp)
    c_2 = JuMP.@NLconstraint(wm.model, r * head_loss(qp) >= inv(L) * dhp)

    # Add constraints for negative flow and head difference.
    qn, dhn = var(wm, n, :qn_pipe, a), var(wm, n, :dhn_pipe, a)
    c_3 = JuMP.@NLconstraint(wm.model, r * head_loss(qn) <= inv(L) * dhn)
    c_4 = JuMP.@NLconstraint(wm.model, r * head_loss(qn) >= inv(L) * dhn)

    # Append the :pipe_head_loss constraint array.
    append!(con(wm, n, :pipe_head_loss)[a], [c_1, c_2, c_3, c_4])
end


"Pump head gain constraint when the pump status is ambiguous."
function constraint_on_off_pump_head_gain(wm::AbstractNCDModel, n::Int, a::Int, node_fr::Int, node_to::Int, pc::Array{Float64}, q_min_forward::Float64)
    # Gather pump flow, head gain, and status variables.
    qp, g, z = var(wm, n, :qp_pump, a), var(wm, n, :g_pump, a), var(wm, n, :z_pump, a)

    # Define the (relaxed) head gain relationship for the pump.
    c_1 = JuMP.@constraint(wm.model, g >= pc[1] * qp^2 + pc[2] * qp + pc[3] * z)
    c_2 = JuMP.@constraint(wm.model, g <= pc[1] * qp^2 + pc[2] * qp + pc[3] * z)
    append!(con(wm, n, :on_off_pump_head_gain, a), [c_1, c_2])
end


"Defines the objective for the owf problem is `NCD` formulations."
function objective_owf_default(wm::AbstractNCDModel)
    objective = zero(JuMP.QuadExpr)

    for (n, nw_ref) in nws(wm)
        efficiency = 0.85 # TODO: How can the efficiency curve be used?
        coeff = _DENSITY * _GRAVITY * ref(wm, n, :time_step) * inv(efficiency)

        for (a, pump) in nw_ref[:pump]
            if haskey(pump, "energy_price")
                # Get pump flow and head gain variables.
                qp, g = var(wm, n, :qp_pump, a), var(wm, n, :g_pump, a)

                # Constrain cost_var and append to the objective expression.
                JuMP.add_to_expression!(objective, coeff * pump["energy_price"], qp, g)
            else
                Memento.error(_LOGGER, "No cost given for pump \"$(pump["name"])\"")
            end
        end
    end

    # Minimize the cost (in units of currency) required to operate pumps.
    return JuMP.@objective(wm.model, _MOI.MIN_SENSE, objective)
end