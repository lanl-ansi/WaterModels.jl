"WaterModels wrapper for the InfrastructureModels `apply!` function."
function apply_wm!(
    func!::Function,
    ref::Dict{Symbol,<:Any},
    data::Dict{String,<:Any};
    apply_to_subnetworks::Bool = true,
)
    _IM.apply!(func!, ref, data, wm_it_sym; apply_to_subnetworks = apply_to_subnetworks)
end


"Filter only the non-inactive components from a dictionary of components."
function _filter_active_components(components::Dict{Int,<:Any})::Dict{Int,<:Any}
    return filter(x -> x.second["status"] !== STATUS_INACTIVE, components)
end


"Build dictionaries that map node indices from and to node-connecting components."
function _build_node_map(
    nodes::Dict{Int,<:Any},
    components::Dict{Int,<:Any},
)::Tuple{Dict,Dict}
    ref_fr = Dict{Int,Array{Int,1}}(i => Array{Int,1}([]) for i in keys(nodes))
    ref_to = Dict{Int,Array{Int,1}}(i => Array{Int,1}([]) for i in keys(nodes))

    for (i, component) in components
        push!(ref_fr[component["node_fr"]], i)
        push!(ref_to[component["node_to"]], i)
    end

    return ref_fr, ref_to
end


"Store pump head gain functional properties in each pump object."
function _set_pump_head_gain_properties!(pumps::Dict{Int,<:Any})
    # Store head curve functional information within each pump object.
    map(x -> x["head_curve_function"] = _calc_head_curve_function(x), values(pumps))
    map(x -> x["head_curve_derivative"] = _calc_head_curve_derivative(x), values(pumps))
    map(x -> x["head_curve_coefficients"] = _calc_head_curve_coefficients(x), values(pumps))
end