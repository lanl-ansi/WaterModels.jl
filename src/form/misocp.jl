# Define MISOCP implementations of Water Models

export
  MISOCPWaterModel, StandardMISOCPForm

""
@compat abstract type AbstractMISOCPForm <: AbstractWaterFormulation end

""
@compat abstract type StandardMISOCPForm <: AbstractMISOCPForm end

const MISOCPWaterModel = GenericWaterModel{StandardMISOCPForm}

"default MISOCP constructor"
MISOCPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMISOCPForm; kwargs...)

""
function variable_flux{T <: AbstractMISOCPForm}(wm::GenericWaterModel{T})
    max_flow = wm.ref[:max_flow]

    wm.var[:l] = @variable(wm.model, [i in keys(wm.ref[:pipe])], basename="l", lowerbound=0.0, upperbound=wm.ref[:connection][i]["resistance"] * max_flow^2, start = getstart(wm.ref[:connection], i, "l_start", 0))
    wm.var[:f] = @variable(wm.model, [i in keys(wm.ref[:connection])], basename="f", lowerbound=-max_flow, upperbound=max_flow, start = getstart(wm.ref[:connection], i, "f_start", 0))

end

# ""
# function variable_flux_ne{T <: AbstractMISOCPForm}(gm::GenericGasModel{T})
#     max_flow = gm.ref[:max_flow]
#     gm.var[:l_ne] = @variable(gm.model, [i in keys(gm.ref[:ne_pipe])], basename="l_ne", lowerbound=0.0, upperbound=1/gm.ref[:ne_connection][i]["resistance"] * max_flow^2, start = getstart(gm.ref[:ne_connection], i, "l_start", 0))
#     gm.var[:f_ne] = @variable(gm.model, [i in keys(gm.ref[:ne_connection])], basename="f_ne", lowerbound=-max_flow, upperbound=max_flow, start = getstart(gm.ref[:ne_connection], i, "f_start", 0))
# end

" Hazen-Williams equation with discrete direction variables "
function constraint_hazen_williams{T <: AbstractMISOCPForm}(wm::GenericWaterModel{T}, pipe_idx)
    pipe = wm.ref[:connection][pipe_idx]
    i_junction_idx = pipe["f_junction"]
    j_junction_idx = pipe["t_junction"]

    hi = wm.var[:h][i_junction_idx]
    hj = wm.var[:h][j_junction_idx]
    yp = wm.var[:yp][pipe_idx]
    yn = wm.var[:yn][pipe_idx]
    l  = wm.var[:l][pipe_idx]
    f  = wm.var[:f][pipe_idx]

    hd_max = pipe["hd_max"]
    hd_min = pipe["hd_min"]
    max_flow = wm.ref[:max_flow]
    # resistance = 4.727*pipe["roughness"]^(-1.852)*(pipe["diameter"])^(-4.871)*pipe["length"];
    c1 = @constraint(wm.model, l >= hj - hi + hd_min*(yp - yn + 1))
    c2 = @constraint(wm.model, l >= hi - hj + hd_max*(yp - yn - 1))
    c3 = @constraint(wm.model, l <= hj - hi + hd_max*(yp - yn + 1))
    c4 = @constraint(wm.model, l <= hi - hj + hd_min*(yp - yn - 1))
    c5 = @constraint(wm.model, l >= pipe["resistance"]*f^2)

    if !haskey(wm.constraint, :hazen_williams1)
        wm.constraint[:hazen_williams1] = Dict{Int,ConstraintRef}()
        wm.constraint[:hazen_williams2] = Dict{Int,ConstraintRef}()
        wm.constraint[:hazen_williams3] = Dict{Int,ConstraintRef}()
        wm.constraint[:hazen_williams4] = Dict{Int,ConstraintRef}()
        wm.constraint[:hazen_williams5] = Dict{Int,ConstraintRef}()
    end
    wm.constraint[:hazen_williams1][pipe_idx] = c1
    wm.constraint[:hazen_williams2][pipe_idx] = c2
    wm.constraint[:hazen_williams3][pipe_idx] = c3
    wm.constraint[:hazen_williams4][pipe_idx] = c4
    wm.constraint[:hazen_williams5][pipe_idx] = c5
end

# " Weisbach equation with discrete direction variables "
# function constraint_weisbach{T <: AbstractMISOCPForm}(wm::GenericWaterModel{T}, pipe_idx)
#     pipe = wm.ref[:connection][pipe_idx]
#     i_junction_idx = pipe["f_junction"]
#     j_junction_idx = pipe["t_junction"]
#
#     hi = wm.var[:h][i_junction_idx]
#     hj = wm.var[:h][j_junction_idx]
#     yp = wm.var[:yp][pipe_idx]
#     yn = wm.var[:yn][pipe_idx]
#     l  = wm.var[:l][pipe_idx]
#     f  = wm.var[:f][pipe_idx]
#
#     hd_max = pipe["hd_max"]
#     hd_min = pipe["hd_min"]
#     max_flow = wm.ref[:max_flow]
#
#     c1 = @constraint(wm.model, l >= hj - hi + hd_min*(yp - yn + 1))
#     c2 = @constraint(wm.model, l >= hi - hj + hd_max*(yp - yn - 1))
#     c3 = @constraint(wm.model, l <= hj - hi + hd_max*(yp - yn + 1))
#     c4 = @constraint(wm.model, l <= hi - hj + hd_min*(yp - yn - 1))
#     c5 = @constraint(wm.model, pipe["resistance"]*l >= f^2)
#
#     if !haskey(wm.constraint, :weisbach1)
#         wm.constraint[:weisbach1] = Dict{Int,ConstraintRef}()
#         wm.constraint[:weisbach2] = Dict{Int,ConstraintRef}()
#         wm.constraint[:weisbach3] = Dict{Int,ConstraintRef}()
#         wm.constraint[:weisbach4] = Dict{Int,ConstraintRef}()
#         wm.constraint[:weisbach5] = Dict{Int,ConstraintRef}()
#     end
#     wm.constraint[:weisbach1][pipe_idx] = c1
#     wm.constraint[:weisbach2][pipe_idx] = c2
#     wm.constraint[:weisbach3][pipe_idx] = c3
#     wm.constraint[:weisbach4][pipe_idx] = c4
#     wm.constraint[:weisbach5][pipe_idx] = c5
# end


# "Weymouth equation with fixed direction"
# function constraint_weymouth_fixed_direction{T <: AbstractMISOCPForm}(gm::GenericGasModel{T}, pipe_idx)
#     pipe = gm.ref[:connection][pipe_idx]
#     i_junction_idx = pipe["f_junction"]
#     j_junction_idx = pipe["t_junction"]
#
#     pi = gm.var[:p][i_junction_idx]
#     pj = gm.var[:p][j_junction_idx]
#     yp = pipe["yp"]
#     yn = pipe["yn"]
#     l  = gm.var[:l][pipe_idx]
#     f  = gm.var[:f][pipe_idx]
#
#     pd_max = pipe["pd_max"]
#     pd_min = pipe["pd_min"]
#     max_flow = gm.ref[:max_flow]
#
#     c1 = @constraint(gm.model, l >= pj - pi + pd_min*(yp - yn + 1))
#     c2 = @constraint(gm.model, l >= pi - pj + pd_max*(yp - yn - 1))
#     c3 = @constraint(gm.model, l <= pj - pi + pd_max*(yp - yn + 1))
#     c4 = @constraint(gm.model, l <= pi - pj + pd_min*(yp - yn - 1))
#     c5 = @constraint(gm.model, pipe["resistance"]*l >= f^2)
#
#     if !haskey(gm.constraint, :weymouth_fixed_direction1)
#         gm.constraint[:weymouth_fixed_direction1] = Dict{Int,ConstraintRef}()
#         gm.constraint[:weymouth_fixed_direction2] = Dict{Int,ConstraintRef}()
#         gm.constraint[:weymouth_fixed_direction3] = Dict{Int,ConstraintRef}()
#         gm.constraint[:weymouth_fixed_direction4] = Dict{Int,ConstraintRef}()
#         gm.constraint[:weymouth_fixed_direction5] = Dict{Int,ConstraintRef}()
#     end
#     gm.constraint[:weymouth_fixed_direction1][pipe_idx] = c1
#     gm.constraint[:weymouth_fixed_direction2][pipe_idx] = c2
#     gm.constraint[:weymouth_fixed_direction3][pipe_idx] = c3
#     gm.constraint[:weymouth_fixed_direction4][pipe_idx] = c4
#     gm.constraint[:weymouth_fixed_direction5][pipe_idx] = c5
# end
