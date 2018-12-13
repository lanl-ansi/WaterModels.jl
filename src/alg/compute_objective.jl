function compute_objective(wm::GenericWaterModel,
                           resistance_indices::Dict{Int, Int},
                           n::Int = wm.cnw)
    objective = 0.0

    for (a, connection) in wm.ref[:nw][n][:connection]
        r_a = resistance_indices[a]
        C_a_r = wm.ref[:nw][wm.cnw][:resistance_cost][a][r_a]
        objective += C_a_r * connection["length"]
    end

    return objective
end
