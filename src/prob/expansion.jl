export run_expansion

function run_expansion(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_expansion; kwargs...)
end

function post_expansion(wm::GenericWaterModel)
    variable_flow(wm)
    variable_head(wm)
    variable_head_difference(wm)
    variable_flow_direction(wm)

    for i in [collect(ids(wm, :junctions)); collect(ids(wm, :reservoirs))]
        constraint_flow_conservation(wm, i)
        constraint_potential_flow_coupling(wm, i)
    end

    for a in collect(ids(wm, :pipes))
        constraint_define_gamma(wm, a)
        constraint_bidirectional_flow(wm, a)
    end
end
