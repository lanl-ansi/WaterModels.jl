function _relax_binary_variable!(v::JuMP.VariableRef)
    JuMP.unset_binary(v) # Make the binary variable continuous and free.

    if !JuMP.is_fixed(v) # If the variable is not fixed, set lower and upper bounds.
        JuMP.set_lower_bound(v, 0.0) # Lower-bound the relaxed binary variables.
        JuMP.set_upper_bound(v, 1.0) # Upper-bound the relaxed binary variables.
    end
end


function relax_all_binary_variables!(wm::AbstractWaterModel)
    vars = filter(v -> JuMP.is_binary(v), JuMP.all_variables(wm.model))
    _relax_binary_variable!.(vars) # Relax all binary variables.
end


function _relax_variables_with_symbol!(wm::AbstractWaterModel, symbol::Symbol)
    for nw in nw_ids(wm) # Loop over all multinetwork subnetworks.
        vars = filter(v -> JuMP.is_binary(v), vcat(var(wm, nw, symbol)...))
        _relax_binary_variable!.(vars)
    end
end


function _relax_all_directions!(wm::AbstractWaterModel)
    _relax_pipe_directions!(wm, :y_pipe)
    _relax_pump_directions!(wm, :y_pump)
    _relax_regulator_directions!(wm, :y_regulator)
    _relax_short_pipe_directions!(wm, :y_short_pipe)
    _relax_valve_directions!(wm, :y_valve)
end


function _relax_all_indicators!(wm::AbstractWaterModel)
    _relax_pump_indicators!(wm, :z_pump)
    _relax_regulator_indicators!(wm, :z_regulator)
    _relax_valve_indicators!(wm, :z_valve)
end
