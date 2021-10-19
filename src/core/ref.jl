"WaterModels wrapper for the InfrastructureModels `apply!` function."
function apply_wm!(func!::Function, ref::Dict{Symbol, <:Any}, data::Dict{String, <:Any}; apply_to_subnetworks::Bool = true)
    _IM.apply!(func!, ref, data, wm_it_sym; apply_to_subnetworks = apply_to_subnetworks)
end