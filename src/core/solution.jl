function _IM.solution_preprocessor(wm::AbstractWaterModel, solution::Dict)
    wm_data = get_wm_data(wm.data)
    solution["it"][wm_it_name]["per_unit"] =
        get_data_wm((x -> return x["per_unit"]),
        wm_data; apply_to_subnetworks = false)

    solution["it"][wm_it_name]["multinetwork"] = ismultinetwork(wm)
    solution["it"][wm_it_name]["base_flow"] = wm.ref[:it][wm_it_sym][:base_flow]
    solution["it"][wm_it_name]["base_head"] = wm.ref[:it][wm_it_sym][:base_head]
    solution["it"][wm_it_name]["base_length"] = wm.ref[:it][wm_it_sym][:base_length]
    solution["it"][wm_it_name]["base_mass"] = wm.ref[:it][wm_it_sym][:base_mass]
    solution["it"][wm_it_name]["base_time"] = wm.ref[:it][wm_it_sym][:base_time]
    solution["it"][wm_it_name]["head_loss"] = wm.ref[:it][wm_it_sym][:head_loss]
end


"WaterModels wrapper for the InfrastructureModels `sol_component_value` function."
function sol_component_value(aim::AbstractWaterModel, n::Int, comp_name::Symbol, field_name::Symbol, comp_ids, variables)
    return _IM.sol_component_value(aim, wm_it_sym, n, comp_name, field_name, comp_ids, variables)
end