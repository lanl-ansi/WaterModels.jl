# Define MINLP implementations of Water Models

export
    MINLPWaterModel, StandardMINLPForm
""
abstract type AbstractMINLPForm <: AbstractWaterFormulation end

""
abstract type StandardMINLPForm <: AbstractMINLPForm end

const MINLPWaterModel = GenericWaterModel{StandardMINLPForm}

"default MINLP constructor"
MINLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMINLPForm; kwargs...)

"Weisbach equation with discrete direction variables "
function constraint_weisbach{T <: AbstractMINLPForm}(wm::GenericWaterModel{T}, pipe_idx)
    pipe = wm.ref[:connection][pipe_idx]
    i_junction_idx = pipe["f_junction"]
    j_junction_idx = pipe["t_junction"]

    i = wm.ref[:junction][i_junction_idx]
    j = wm.ref[:junction][j_junction_idx]

    hi = wm.var[:h][i_junction_idx]
    hj = wm.var[:h][j_junction_idx]
    yp = wm.var[:yp][pipe_idx]
    yn = wm.var[:yn][pipe_idx]
    f  = wm.var[:f][pipe_idx]

    max_flow = wm.ref[:max_flow]
    w = pipe["resistance"]

    c1 = @NLconstraint(wm.model, w*(hi - hj) >= f^2 - (1-yp)*max_flow^2)
    c2 = @NLconstraint(wm.model, w*(hi - hj) <= f^2 + (1-yp)*max_flow^2)
    c3 = @NLconstraint(wm.model, w*(hj - hi) >= f^2 - (1-yn)*max_flow^2)
    c4 = @NLconstraint(wm.model, w*(hj - hi) <= f^2 + (1-yn)*max_flow^2)

   if !haskey(wm.constraint, :weisbach1)
        wm.constraint[:weisbach1] = Dict{Int,ConstraintRef}()
        wm.constraint[:weisbach2] = Dict{Int,ConstraintRef}()
        wm.constraint[:weisbach3] = Dict{Int,ConstraintRef}()
        wm.constraint[:weisbach4] = Dict{Int,ConstraintRef}()
    end
    wm.constraint[:weisbach1][pipe_idx] = c1
    wm.constraint[:weisbach2][pipe_idx] = c2
    wm.constraint[:weisbach3][pipe_idx] = c3
    wm.constraint[:weisbach4][pipe_idx] = c4
end

# "Weymouth equation with fixed direction variables"
# function constraint_weymouth_fixed_direction{T <: AbstractMINLPForm}(gm::GenericGasModel{T}, pipe_idx)
#     pipe = gm.ref[:connection][pipe_idx]
#     i_junction_idx = pipe["f_junction"]
#     j_junction_idx = pipe["t_junction"]
#
#     i = gm.ref[:junction][i_junction_idx]
#     j = gm.ref[:junction][j_junction_idx]
#
#     pi = gm.var[:p][i_junction_idx]
#     pj = gm.var[:p][j_junction_idx]
#     yp = pipe["yp"]
#     yn = pipe["yn"]
#     f  = gm.var[:f][pipe_idx]
#
#     max_flow = gm.ref[:max_flow]
#     w = pipe["resistance"]
#
#     c1 = @NLconstraint(gm.model, w*(pi - pj) >= f^2 - (1-yp)*max_flow^2)
#     c2 = @NLconstraint(gm.model, w*(pi - pj) <= f^2 + (1-yp)*max_flow^2)
#     c3 = @NLconstraint(gm.model, w*(pj - pi) >= f^2 - (1-yn)*max_flow^2)
#     c4 = @NLconstraint(gm.model, w*(pj - pi) <= f^2 + (1-yn)*max_flow^2)
#
#    if !haskey(gm.constraint, :weymouth_fixed_direction1)
#         gm.constraint[:weymouth_fixed_direction1] = Dict{Int,ConstraintRef}()
#         gm.constraint[:weymouth_fixed_direction2] = Dict{Int,ConstraintRef}()
#         gm.constraint[:weymouth_fixed_direction3] = Dict{Int,ConstraintRef}()
#         gm.constraint[:weymouth_fixed_direction4] = Dict{Int,ConstraintRef}()
#     end
#     gm.constraint[:weymouth_fixed_direction1][pipe_idx] = c1
#     gm.constraint[:weymouth_fixed_direction2][pipe_idx] = c2
#     gm.constraint[:weymouth_fixed_direction3][pipe_idx] = c3
#     gm.constraint[:weymouth_fixed_direction4][pipe_idx] = c4
# end
