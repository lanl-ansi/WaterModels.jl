# Defines MINLP implementations of water distribution models.

export MINLPWaterModel, StandardMINLPForm

abstract type AbstractMINLPForm <: AbstractWaterFormulation end
abstract type StandardMINLPForm <: AbstractMINLPForm end

"Default (nonconvex) MINLP model."
const MINLPWaterModel = GenericWaterModel{StandardMINLPForm}

"Default MINLP constructor."
MINLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMINLPForm; kwargs...)

function variable_flow(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw; alpha::Float64=1.852) where T <: AbstractMINLPForm
    variable_directed_flow(wm, n, alpha=alpha, bounded=true)
    variable_directed_flow_ne(wm, n, alpha=alpha, bounded=true)
end

function constraint_potential_loss(wm::GenericWaterModel{T}, a::Int, n::Int) where T <: StandardMINLPForm
    if !haskey(wm.con[:nw][n], :potential_loss_1)
        wm.con[:nw][n][:potential_loss_1] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        wm.con[:nw][n][:potential_loss_2] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    wm.con[:nw][n][:potential_loss_1][a] = Dict{Int, JuMP.ConstraintRef}()
    wm.con[:nw][n][:potential_loss_2][a] = Dict{Int, JuMP.ConstraintRef}()

    dhp = wm.var[:nw][n][:dhp][a]
    dhn = wm.var[:nw][n][:dhn][a]
    L = wm.ref[:nw][n][:links][a]["length"]

    for r in 1:length(wm.ref[:nw][n][:resistance][a])
        q_n_a_r = wm.var[:nw][n][:qn][a][r]
        q_p_a_r = wm.var[:nw][n][:qp][a][r]
        resistance = wm.ref[:nw][n][:resistance][a][r]

        con_1 = JuMP.@NLconstraint(wm.model, dhp / L >= resistance * f_alpha(q_p_a_r))
        wm.con[:nw][n][:potential_loss_1][a][r] = con_1

        con_2 = JuMP.@NLconstraint(wm.model, dhn / L >= resistance * f_alpha(q_n_a_r))
        wm.con[:nw][n][:potential_loss_2][a][r] = con_2
    end
end


#function constraint_flow_direction(wm::GenericWaterModel{T}, a::Int, n::Int = wm.cnw) where T <: StandardMINLPForm
#    # Get source and target nodes corresponding to the edge.
#    i = parse(Int, wm.ref[:nw][n][:pipes][a]["node1"])
#    j = parse(Int, wm.ref[:nw][n][:pipes][a]["node2"])
#
#    # Collect variables needed for the constraint.
#    q = wm.var[:nw][n][:q][a]
#    h_i = wm.var[:nw][n][:h][i]
#    h_j = wm.var[:nw][n][:h][j]
#    dhp = wm.var[:nw][n][:dhp][a]
#    dhn = wm.var[:nw][n][:dhn][a]
#end

#"Non-convex Darcy-Weisbach constraint with unknown direction."
#function constraint_dw_unknown_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMINLPForm
#    # Collect variables and parameters needed for the constraint.
#    q, h_i, h_j, viscosity, lambda = get_dw_requirements(wm, a, n)
#
#    # Add a nonlinear constraint for the head loss.
#    @NLconstraint(wm.model, h_i - h_j == lambda * q * abs(q))
#end
#
#"Non-convex Hazen-Williams constraint for flow with unknown direction."
#function constraint_hw_unknown_direction(wm::GenericWaterModel{T}, a, n::Int = wm.cnw) where T <: StandardMINLPForm
#    # Collect variables and parameters needed for the constraint.
#    q, h_i, h_j, lambda = get_hw_requirements(wm, a, n)
#
#    # Add a non-convex constraint for the head loss.
#    @NLconstraint(wm.model, h_i - h_j == lambda * q * (q^2)^0.426)
#end
#
