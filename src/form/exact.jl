"These problem forms use exact versions of the head loss equation when the flow direction is fixed."
AbstractExactForm = Union{AbstractMINLPBForm, AbstractNLPForm}

"Exact Hazen-Williams constraint for flow with known direction."
function constraint_hw_known_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: AbstractExactForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a non-convex constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) == lambda * (q^2)^0.926)
end

"Exact Darcy-Weisbach constraint with known direction."
function constraint_dw_known_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: AbstractExactForm
    # Collect variables and parameters needed for the constraint.
    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)

    # Fix the direction associated with the flow.
    dir = Int(wm.data["pipes"][a]["flow_direction"])
    fix_flow_direction(q, dir)

    # Add a non-convex quadratic constraint for the head loss.
    @NLconstraint(wm.model, dir * (h_i - h_j) == lambda * q^2)
end
