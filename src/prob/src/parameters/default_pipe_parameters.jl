# default_network_parameters
default_pipe_data = Dict()

get!(default_pipe_data, "steam_pipe_friction_factor", 0.01)
get!(default_pipe_data, "water_pipe_friction_factor", 0.01)
get!(default_pipe_data, "steam_pipe_thermal_loss_coefficient", 10)
get!(default_pipe_data, "water_pipe_thermal_loss_coefficient", 1e-5)
