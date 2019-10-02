# Define MICP (mixed-integer convex program) implementations of water distribution models.

function variable_pump(wm::AbstractMICPModel, n::Int=wm.cnw) end
function constraint_head_loss_pump(wm::AbstractMICPModel, n::Int, a::Int, node_fr::Int, node_to::Int) end

function constraint_head_loss_pipe_ne(wm::AbstractMICPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, pipe_resistances) 
    for (r_id, r) in enumerate(pipe_resistances)
        dhn = var(wm, n, :dhn, a)
        qn_ne = var(wm, n, :qn_ne, a)[r_id]
        cn = JuMP.@NLconstraint(wm.model, r * head_loss(qn_ne) <= inv(L) * dhn)

        dhp = var(wm, n, :dhp, a)
        qp_ne = var(wm, n, :qp_ne, a)[r_id]
        cp = JuMP.@NLconstraint(wm.model, r * head_loss(qp_ne) <= inv(L) * dhp)

        append!(con(wm, n, :head_loss, a), [cn, cp])
    end
end

function constraint_head_loss_pipe(wm::AbstractMICPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    qp = var(wm, n, :qp, a)
    dhp = var(wm, n, :dhp, a)
    cp = JuMP.@NLconstraint(wm.model, r * head_loss(qp) <= inv(L) * dhp)

    qn = var(wm, n, :qn, a)
    dhn = var(wm, n, :dhn, a)
    cn = JuMP.@NLconstraint(wm.model, r * head_loss(qn) <= inv(L) * dhn)

    append!(con(wm, n, :head_loss)[a], [cp, cn])
end

#function objective_wf(wm::AbstractMICPModel)
#    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
#end
#
#function objective_owf(wm::AbstractMICPModel) 
#    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
#end
