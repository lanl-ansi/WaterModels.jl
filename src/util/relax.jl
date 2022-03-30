function _relax_binary_variable!(v::JuMP.VariableRef)
    if JuMP.is_binary(v)
        JuMP.unset_binary(v) # Make the binary variable continuous and free.
    end

    if !JuMP.is_fixed(v) # If the variable is not fixed, set lower and upper bounds.
        JuMP.set_lower_bound(v, 0.0) # Lower-bound the relaxed binary variables.
        JuMP.set_upper_bound(v, 1.0) # Upper-bound the relaxed binary variables.
    end
end


function set_binary_variables!(vars::Array{JuMP.VariableRef, 1})
    map(x -> JuMP.set_binary(x), vars)
end


function get_all_binary_vars_at_nw!(wm::AbstractWaterModel, nw::Int)
    vars_binary = Array{JuMP.VariableRef, 1}([])
    
    for var_entry in values(var(wm, nw))
        vars = filter(x -> isa(x, JuMP.VariableRef), vcat(var_entry...))
        vars_binary_inner = collect(filter(x -> JuMP.is_binary(x), vars))
        append!(vars_binary, vars_binary_inner)
    end

    return vars_binary
end



function relax_all_binary_variables_at_nw!(wm::AbstractWaterModel, nw::Int)
    vars = get_all_binary_vars_at_nw!(wm, nw)
    map(x -> JuMP.unset_binary(x), vars)
end


function relax_all_binary_variables!(wm::AbstractWaterModel)
    vars = filter(v -> JuMP.is_binary(v), JuMP.all_variables(wm.model))
    _relax_binary_variable!.(vars) # Relax all binary variables.
end


function _relax_variables_with_symbol!(wm::AbstractWaterModel, symbol::Symbol)
    for nw in sort(collect(nw_ids(wm)))[1:end-1]
        vars = filter(v -> JuMP.is_binary(v), vcat(var(wm, nw, symbol)...))
        _relax_binary_variable!.(vars)
    end
end


function _relax_all_direction_variables!(wm::AbstractWaterModel)
    _relax_variables_with_symbol!(wm, :y_pipe)
    _relax_variables_with_symbol!(wm, :y_pump)
    _relax_variables_with_symbol!(wm, :y_regulator)
    _relax_variables_with_symbol!(wm, :y_short_pipe)
    _relax_variables_with_symbol!(wm, :y_ne_short_pipe)
    _relax_variables_with_symbol!(wm, :y_valve)
end


function _relax_all_indicator_variables!(wm::AbstractWaterModel)
    _relax_variables_with_symbol!(wm, :z_pump)
    _relax_variables_with_symbol!(wm, :z_regulator)
    _relax_variables_with_symbol!(wm, :z_valve)
end


function relax_every_other_indicator_variable!(wm::AbstractWaterModel, step::Int)
    var_symbols = Array{Symbol}([:z_pump, :z_regulator, :z_valve])

    for nw in sort(collect(nw_ids(wm)))[2:step:end-1]
        vars = vcat([vcat(var(wm, nw, s)...) for s in var_symbols]...)
        _relax_binary_variable!.(vars)
    end
end


function _relax_last_indicator_variables!(wm::AbstractWaterModel; last_num_steps::Int = length(nw_ids(wm)))
    var_symbols = Array{Symbol}([:z_pump, :z_regulator, :z_valve])
    network_ids = reverse(sort(collect(nw_ids(wm)))[1:end-1])

    for nw in network_ids[1:min(length(network_ids), last_num_steps)]
        vars = vcat([vcat(var(wm, nw, s)...) for s in var_symbols]...)
        _relax_binary_variable!.(vars)
    end
end


function _relax_last_variables!(wm::AbstractWaterModel; last_num_steps::Int = length(nw_ids(wm)))
    var_symbols = Array{Symbol}([:z_pump, :z_regulator, :z_valve])
    network_ids = reverse(sort(collect(nw_ids(wm)))[1:end-1])

    for nw in network_ids[1:min(length(network_ids), last_num_steps)]
        vars = vcat([vcat(var(wm, nw, s)...) for s in var_symbols]...)
        _relax_binary_variable!.(vars)
    end
end