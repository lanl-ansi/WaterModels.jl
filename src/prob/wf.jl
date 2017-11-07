# Definitions for running a feasible water flow

# Note that this particular formulation assumes the binary variable implementation of flow direction
# We would need to do some abstraction to support the absolute value formulation

export run_wf

# entry point into running the water flow feasability problem
function run_wf(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_wf; kwargs...)
end

""
function run_soc_wf(file, solver; kwargs...)
    return run_wf(file, MISOCPWaterModel, solver; kwargs...)
end

# construct the water flow feasbility problem
function post_wf{T}(wm::GenericWaterModel{T})
    # variable_pressure_sqr(wm)

    variable_head(wm)
    # variable_head(wm)
    variable_flux(wm)
    variable_connection_direction(wm)
    # variable_flux_square(wm)
    variable_valve_operation(wm)


    for (i,junction) in wm.ref[:junction]

        if junction["type"] == "regular"
          constraint_junction_flow_balance(wm, i)
          constraint_elevation_bound_regular(wm, i)
        end

        # if junction["qgfirm"] > 0.0 && junction["qlfirm"] == 0.0
        if junction["type"] == "reservoir"
            constraint_source_flow(wm, i)
            constraint_elevation_bound_source(wm, i)
        end

        # if junction["qgfirm"] == 0.0 && junction["qlfirm"] > 0.0
        #     constraint_sink_flow(wm, i)
        # end

        if junction["type"] == "regular" && junction["demand"] == 0.0 && junction["degree"] == 2
           constraint_conserve_flow(wm, i)
        end
    end

    for (i,connection) in wm.ref[:connection]
        constraint_flow_direction_choice(wm, i)
        # constraint_parallel_flow(wm, i)
    end

    # for i in [collect(keys(wm.ref[:pipe])); collect(keys(wm.ref[:resistor]))]
    for i in keys(wm.ref[:pipe])
        constraint_on_off_head_drop(wm, i)
        constraint_on_off_pipe_flow_direction(wm, i)
        # constraint_on_off_pipe_flow_direction_diameter(wm, i)
        constraint_hazen_williams(wm, i)
    end

    # for (i,pipe) in wm.ref[:short_pipe]
    #     constraint_short_pipe_pressure_drop(wm, i)
    #     constraint_on_off_short_pipe_flow_direction(wm, i)
    # end

    # for (i, compressor) in wm.ref[:compressor]
    #     constraint_on_off_compressor_flow_direction(wm, i)
    #     constraint_on_off_compressor_ratios(wm, i)
    # end

    for (i,valve) in wm.ref[:valve]
        constraint_on_off_valve_flow_direction(wm, i)
        constraint_on_off_valve_pressure_drop(wm, i)
    end

    for (i, valve) in wm.ref[:control_valve]
        constraint_on_off_control_valve_flow_direction(wm, i)
        constraint_on_off_control_valve_pressure_drop(wm, i)
    end
end
