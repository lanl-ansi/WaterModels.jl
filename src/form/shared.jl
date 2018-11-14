"Set new bounds for q given some specified direction of flow (-1 or 1)."
function fix_flow_direction(q::JuMP.Variable, direction::Int)
    # Fix the direction of the flow.
    setlowerbound(q, direction == 1 ? 0.0 : getlowerbound(q))
    setupperbound(q, direction == 1 ? getupperbound(q) : 0.0)
end

"Get variables commonly used in the construction of head loss constraints."
function get_common_variables(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: AbstractWaterFormulation
    # Get source and target nodes corresponding to the edge.
    i = parse(Int, wm.ref[:nw][n][:pipes][a]["node1"])
    j = parse(Int, wm.ref[:nw][n][:pipes][a]["node2"])

    # Collect variables needed for the constraint.
    q = wm.var[:nw][n][:q][a]
    h_i = wm.var[:nw][n][:h][i]
    h_j = wm.var[:nw][n][:h][j]

    # Return the variables.
    return q, h_i, h_j
end

"Get variables and constants used in the construction of Darcy-Weisbach constraints."
function get_dw_requirements(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: AbstractWaterFormulation
    q, h_i, h_j = get_common_variables(wm, a, n)
    viscosity = wm.ref[:nw][n][:options]["viscosity"]
    lambda = calc_friction_factor_dw(wm.ref[:nw][n][:pipes][a], viscosity)
    return q, h_i, h_j, viscosity, lambda
end

"Get variables and constants used in the construction of Hazen-Williams constraints."
function get_hw_requirements(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: AbstractWaterFormulation
    q, h_i, h_j = get_common_variables(wm, a, n)
    lambda = calc_friction_factor_hw(wm.ref[:nw][n][:pipes][a])
    return q, h_i, h_j, lambda
end
