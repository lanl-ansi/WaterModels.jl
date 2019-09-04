function relax_integrality!(wm::AbstractWaterModel)
    for var in filter(v -> JuMP.is_binary(v), JuMP.all_variables(wm.model))
        JuMP.unset_binary(var)
        JuMP.set_lower_bound(var, 0.0)
        JuMP.set_upper_bound(var, 1.0)
    end
end
