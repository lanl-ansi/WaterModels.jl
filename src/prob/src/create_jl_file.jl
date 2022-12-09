function create_jl_file(data, default_pipe_data, model_name::AbstractString)
    open(string("network_files/input_",model_name,".jl"), "w") do f
    junction_data = data["junctions"]
    pipes_data = data["pipes"]
    plants_data = data["plants"]
    loads_data = data["loads"]

    lambda_s = default_pipe_data["steam_pipe_friction_factor"]
    lambda_w = default_pipe_data["water_pipe_friction_factor"]
    gamma_s = default_pipe_data["steam_pipe_thermal_loss_coefficient"]
    gamma_w = default_pipe_data["water_pipe_thermal_loss_coefficient"]

    write(f, "junctions_mat = [ \n")
    for junction in junction_data
        write(f, "$(junction["id"]) '$(junction["type"])' $(junction["x_location"]) $(junction["y_location"])\n")
    end
    write(f, "];\n")

    write(f, "pipes_mat = [ \n")
    for pipe in pipes_data
        if(pipe["type"] == 'o')
            lambda = lambda_s
            gamma = gamma_s
        elseif(pipe["type"] == 'r')
            lambda = lambda_w
            gamma = gamma_w
        end
        write(f, "$(pipe["id"]) '$(pipe["type"])' $(pipe["fr_node"]) $(pipe["to_node"]) $(pipe["diameter"]) $(pipe["length"]) $(lambda) $(gamma)\n")
    end
    write(f, "];\n")

    write(f, "plants_mat = [ \n")
    for plant in plants_data
        write(f, "$(plant["id"]) $(plant["fr_node"]) $(plant["to_node"]) $(plant["load"])\n")
    end
    write(f, "];\n")


    write(f, "loads_mat = [ \n")
    for load in loads_data
        write(f, "$(load["id"]) $(load["fr_node"]) $(load["to_node"]) $(load["load"])\n")
    end
    write(f, "];\n")
    close(f)
end

end
