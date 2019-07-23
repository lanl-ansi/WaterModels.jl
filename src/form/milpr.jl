# Define MILPR (mixed-integer linear, relaxed program) implementations of water distribution models.
export MILPRWaterModel, StandardMILPRForm

abstract type AbstractMILPRForm <: AbstractWaterFormulation end
abstract type StandardMILPRForm <: AbstractMILPRForm end

"The default MILPR (mixed-integer linear, relaxed) model is a linear outer-approximation of the MICP model."
const MILPRWaterModel = GenericWaterModel{StandardMILPRForm}

"Default MILPR constructor."
MILPRWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMILPRForm; kwargs...)

function variable_head(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractMILPRForm
    variable_pressure_head(wm, n)
    variable_directed_head_difference(wm, n)
end

function variable_flow(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractMILPRForm
    variable_undirected_flow(wm, n, bounded=true)
    variable_directed_flow(wm, n, bounded=true)
    variable_flow_direction(wm, n)
end

function variable_flow_ne(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractMILPRForm
    variable_directed_flow_ne(wm, n, bounded=true)
end

function variable_pump(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: AbstractMILPRForm
end

function constraint_resistance_selection_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractMILPRForm
    constraint_directed_resistance_selection_ne(wm, a, n)
end

function constraint_potential_loss(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractMILPRForm
    constraint_directed_head_difference(wm, a, n)
    constraint_flow_direction_selection(wm, a, n)
    constraint_directed_potential_loss_ub(wm, a, n)
    constraint_directed_potential_loss(wm, a, n)
end

function constraint_potential_loss_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractMILPRForm
    constraint_directed_head_difference(wm, a, n)
    constraint_flow_direction_selection_ne(wm, a, n)
    constraint_directed_potential_loss_ub_ne(wm, a, n)
    constraint_directed_potential_loss_ne(wm, a, n)
end

function constraint_flow_conservation(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractMILPRForm
    constraint_directed_flow_conservation(wm, i, n)
end

function constraint_flow_conservation_ne(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractMILPRForm
    constraint_directed_flow_conservation_ne(wm, i, n)
end

function constraint_link_flow_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractMILPRForm
    constraint_link_directed_flow_ne(wm, a, n)
end

function constraint_link_flow(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractMILPRForm
    constraint_link_directed_flow(wm, a, n)
end

function constraint_source_flow(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractMILPRForm
    constraint_directed_source_flow(wm, i, n)
end

function constraint_sink_flow(wm::GenericWaterModel{T}, i::Int, n::Int=wm.cnw) where T <: AbstractMILPRForm
    constraint_directed_sink_flow(wm, i, n)
end

function get_linear_outer_approximation(q::JuMP.VariableRef, q_hat::Float64, alpha::Float64)
    return q_hat^alpha + alpha * q_hat^(alpha - 1.0) * (q - q_hat)
end

function constraint_directed_potential_loss_ne(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractMILPRForm
    if !haskey(con(wm, n), :potential_loss_n_ne)
        con(wm, n)[:potential_loss_n_ne] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        con(wm, n)[:potential_loss_p_ne] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    end

    con(wm, n, :potential_loss_n_ne)[a] = Dict{Int, JuMP.ConstraintRef}()
    con(wm, n, :potential_loss_p_ne)[a] = Dict{Int, JuMP.ConstraintRef}()

    alpha = ref(wm, n, :alpha)
    L = ref(wm, n, :links, a)["length"]

    for (r_id, r) in enumerate(ref(wm, n, :resistance, a))
        qn_ne = var(wm, n, :qn_ne, a)[r_id]
        qn_ne_ub = JuMP.upper_bound(qn_ne)
        dhn = var(wm, n, :dhn, a)

        if qn_ne_ub > 0.0
            for q_hat in range(0.0, stop=qn_ne_ub, length=5)
                cut_lhs = r * get_linear_outer_approximation(qn_ne, q_hat, alpha)
                con_n = JuMP.@constraint(wm.model, cut_lhs - inv(L) * dhn <= 0.0)
            end
        end

        qp_ne = var(wm, n, :qp_ne, a)[r_id]
        qp_ne_ub = JuMP.upper_bound(qp_ne)
        dhp = var(wm, n, :dhp, a)

        if qp_ne_ub > 0.0
            for q_hat in range(0.0, stop=qp_ne_ub, length=5)
                cut_lhs = r * get_linear_outer_approximation(qp_ne, q_hat, alpha)
                con_p = JuMP.@constraint(wm.model, cut_lhs - inv(L) * dhp <= 0.0)
            end
        end
    end
end

function constraint_directed_potential_loss(wm::GenericWaterModel{T}, a::Int, n::Int=wm.cnw) where T <: AbstractMILPRForm
    if !haskey(con(wm, n), :potential_loss_n)
        con(wm, n)[:potential_loss_n] = Dict{Int, JuMP.ConstraintRef}()
        con(wm, n)[:potential_loss_p] = Dict{Int, JuMP.ConstraintRef}()
    end

    alpha = ref(wm, n, :alpha)
    L = ref(wm, n, :links, a)["length"]
    r = minimum(ref(wm, n, :resistance, a))

    qn = var(wm, n, :qn, a)
    qn_ub = JuMP.upper_bound(qn)
    dhn = var(wm, n, :dhn, a)

    if qn_ub > 0.0
        for q_hat in range(0.0, stop=qn_ub, length=5)
            cut_lhs = r * get_linear_outer_approximation(qn, q_hat, alpha)
            con_n = JuMP.@constraint(wm.model, cut_lhs - inv(L) * dhn <= 0.0)
        end
    end

    qp = var(wm, n, :qp, a)
    qp_ub = JuMP.upper_bound(qp)
    dhp = var(wm, n, :dhp, a)

    if qp_ub > 0.0
        for q_hat in range(0.0, stop=qp_ub, length=5)
            cut_lhs = r * get_linear_outer_approximation(qp, q_hat, alpha)
            con_p = JuMP.@constraint(wm.model, cut_lhs - inv(L) * dhp <= 0.0)
        end
    end
end

function objective_wf(wm::GenericWaterModel{T}, n::Int=wm.cnw) where T <: StandardMILPRForm
    JuMP.set_objective_sense(wm.model, MOI.FEASIBILITY_SENSE)
end
