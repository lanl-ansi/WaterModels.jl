function unbinarize_mn(wm::AbstractWaterModel, num_steps::Int)
    var_symbols = [:z_check_valve, :z_shutoff_valve, :z_regulator, :z_pump]
    network_ids = sort(collect(nw_ids(wm)), rev=true) # Descending indices.

    for nw in network_ids[1:min(length(network_ids), num_steps)]
        vars = vcat([vcat(var(wm, nw, s)[:]...) for s in var_symbols]...)
        binary_vars = filter(v -> JuMP.is_binary(v), vars)
        JuMP.unset_binary.(binary_vars)
        JuMP.set_lower_bound.(binary_vars, 0.0)
        JuMP.set_upper_bound.(binary_vars, 1.0)
    end
end
