# Constraint templates help simplify data wrangling across multiple Water Flow
# formulations by providing an abstraction layer between the network data and
# network constraint definitions. The constraint template's job is to extract
# the required parameters from a given network data structure and pass the data
# as named arguments to the Water Flow formulations.
#
# Constraint templates should always be defined over "GenericWaterModel" and
# should never refer to model variables.

""
function constraint_tank_state(wm::GenericWaterModel, i::Int; nw::Int=wm.cnw)
    tank = ref(wm, nw, :tanks, i)
    time_step = ref(wm, nw, :options, "time")["hydraulic_timestep"]

    if time_step <= 0.0
        Memento.error(_LOGGER, "Tank states cannot be controlled without a time step.")
    end

    initial_level = ref(wm, nw, :tanks, i)["init_level"]
    surface_area = 0.25 * pi * ref(wm, nw, :tanks, i)["diameter"]^2
    V_initial = surface_area * initial_level

    constraint_tank_state_initial(wm, nw, i, V_initial, convert(Float64, time_step))
end

function constraint_tank_state(wm::GenericWaterModel, i::Int, nw_1::Int, nw_2::Int)
    tank = ref(wm, nw_2, :tanks, i)

    if haskey(ref(wm, nw), :time_series)
        time_step = ref(wm, nw, :time_series)["time_step"]
    else
        Memento.error(_LOGGER, "Tank states cannot be controlled outside of a time series.")
    end

    constraint_tank_state(wm, nw_1, nw_2, i, time_step)
end
