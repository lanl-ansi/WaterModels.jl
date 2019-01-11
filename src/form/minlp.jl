# Defines MINLP implementations of water distribution models.

export MINLPWaterModel, StandardMINLPForm

abstract type AbstractMINLPForm <: AbstractWaterFormulation end
abstract type StandardMINLPForm <: AbstractMINLPForm end

"Default (nonconvex) MINLP model."
const MINLPWaterModel = GenericWaterModel{StandardMINLPForm}

"Default MINLP constructor."
MINLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMINLPForm; kwargs...)

function constraint_potential_loss_segmented(wm::GenericWaterModel{T}, a::Int, n_n::Int, n_s::Int) where T <: StandardMINLPForm
    if !haskey(wm.con[:nw][n_n], :potential_loss_1)
        function_head_loss_hw(wm) # Register the head loss JuMP function.
        wm.con[:nw][n_n][:potential_loss_1] = Dict{Int, Dict{Int, Dict{Int, ConstraintRef}}}()
        wm.con[:nw][n_n][:potential_loss_2] = Dict{Int, Dict{Int, Dict{Int, ConstraintRef}}}()
    end

    wm.con[:nw][n_n][:potential_loss_1][a] = Dict{Int, Dict{Int, ConstraintRef}}()
    wm.con[:nw][n_n][:potential_loss_2][a] = Dict{Int, Dict{Int, ConstraintRef}}()

    dhp = wm.var[:nw][n_n][:dhp][a]
    dhn = wm.var[:nw][n_n][:dhn][a]
    L = wm.ref[:nw][n_n][:connection][a]["length"]

    for k in 1:n_s
        wm.con[:nw][n_n][:potential_loss_1][a][k] = Dict{Int, ConstraintRef}()
        wm.con[:nw][n_n][:potential_loss_2][a][k] = Dict{Int, ConstraintRef}()

        for r in 1:length(wm.ref[:nw][n_n][:resistance][a])
            q_n_akr = wm.var[:nw][n_n][:qn][a][k, r]
            q_p_akr = wm.var[:nw][n_n][:qp][a][k, r]
            resistance = wm.ref[:nw][n_n][:resistance][a][r]

            con_1 = @NLconstraint(wm.model, dhn / L >= resistance * head_loss_hw(q_n_akr))
            wm.con[:nw][n_n][:potential_loss_1][a][k][r] = con_1

            con_2 = @NLconstraint(wm.model, dhp / L >= resistance * head_loss_hw(q_p_akr))
            wm.con[:nw][n_n][:potential_loss_2][a][k][r] = con_2
        end
    end
end

function constraint_potential_loss(wm::GenericWaterModel{T}, a::Int, n_n::Int) where T <: StandardMINLPForm
    if !haskey(wm.con[:nw][n_n], :potential_loss_1)
        function_head_loss_hw(wm) # Register the head loss JuMP function.
        wm.con[:nw][n_n][:potential_loss_1] = Dict{Int, Dict{Int, ConstraintRef}}()
        wm.con[:nw][n_n][:potential_loss_2] = Dict{Int, Dict{Int, ConstraintRef}}()
    end

    wm.con[:nw][n_n][:potential_loss_1][a] = Dict{Int, ConstraintRef}()
    wm.con[:nw][n_n][:potential_loss_2][a] = Dict{Int, ConstraintRef}()

    dhp = wm.var[:nw][n_n][:dhp][a]
    dhn = wm.var[:nw][n_n][:dhn][a]
    L = wm.ref[:nw][n_n][:connection][a]["length"]

    for r in 1:length(wm.ref[:nw][n_n][:resistance][a])
        q_n_r = wm.var[:nw][n_n][:qn][a][r]
        q_p_r = wm.var[:nw][n_n][:qp][a][r]
        resistance = wm.ref[:nw][n_n][:resistance][a][r]

        con_1 = @NLconstraint(wm.model, dhp / L >= resistance * head_loss_hw(q_p_r))
        wm.con[:nw][n_n][:potential_loss_1][a][r] = con_1

        con_2 = @NLconstraint(wm.model, dhn / L >= resistance * head_loss_hw(q_n_r))
        wm.con[:nw][n_n][:potential_loss_2][a][r] = con_2
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
