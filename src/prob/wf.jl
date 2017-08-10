# Definitions for running a feasible water flow

# Note that this particular formulation assumes the binary variable implementation of flow direction
# We would need to do some abstraction to support the absolute value formulation

export run_wf

# entry point into running the water flow feasability problem
function run_wf(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_wf; kwargs...)
end

# construct the water flow feasbility problem
function post_wf{T}(wm::GenericWaterModel{T})
    variable_pressure_sqr(wm)
    variable_flux(wm)
    variable_connection_direction(wm)
    variable_flux_square(wm)
    variable_valve_operation(wm)

    for (i,junction) in wm.set.junctions
        constraint_junction_flow_balance(wm, junction)

        if junction["qgfirm"] > 0.0 && junction["qlfirm"] == 0.0
            constraint_source_flow(wm, junction)
        end

        if junction["qgfirm"] == 0.0 && junction["qlfirm"] > 0.0
            constraint_sink_flow(wm, junction)
        end

        if junction["qgfirm"] == 0.0 && junction["qlfirm"] == 0.0 && junction["degree"] == 2
           constraint_conserve_flow(wm, junction)
        end

    end

    for (i,connection) in wm.set.connections
        constraint_flow_direction_choice(wm, connection)
        constraint_parallel_flow(wm,connection)
    end

    for i in [wm.set.pipe_indexes; wm.set.resistor_indexes]
        pipe = wm.set.connections[i]
        constraint_on_off_pressure_drop(wm, pipe)
        constraint_on_off_pipe_flow_direction(wm,pipe)
        constraint_weisbach(wm,pipe)
    end

    for i in wm.set.short_pipe_indexes
        pipe = wm.set.connections[i]
        constraint_short_pipe_pressure_drop(wm, pipe)
        constraint_on_off_short_pipe_flow_direction(wm,pipe)
    end

    # for i in wm.set.compressor_indexes
    #     compressor = wm.set.connections[i]
    #     constraint_on_off_compressor_flow_direction(wm, compressor)
    #     constraint_on_off_compressor_ratios(wm, compressor)
    # end

    for i in wm.set.valve_indexes
        valve = wm.set.connections[i]
        constraint_on_off_valve_flow_direction(wm, valve)
        constraint_on_off_valve_pressure_drop(wm, valve)
    end

    for i in wm.set.control_valve_indexes
        valve = wm.set.connections[i]
        constraint_on_off_control_valve_flow_direction(wm, valve)
        constraint_on_off_control_valve_pressure_drop(wm, valve)
    end

end
