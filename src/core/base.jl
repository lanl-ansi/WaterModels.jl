"Root of the water formulation type hierarchy."
abstract type AbstractWaterModel <: _IM.AbstractInfrastructureModel end


"Macro for adding base WaterModels fields to the type definition."
_IM.@def wm_fields begin
    WaterModels.@im_fields
end


"""
    solve_model(file::String, model_type::Type, optimizer::Any, build_method::Function; kwargs...)::Dict{String,<:Any}

Instantiates and solves the modeling object from input file with file path `file`. Here,
`model_type` is the model formulation type (e.g., NCWaterModel), `optimizer` is the
optimization solver used to solve the problem (e.g., Gurobi.Optimizer), and
`build_method` is the function used for building the problem specification being
considered (e.g., build_mn_owf). Returns a dictionary of optimization results.
"""
function solve_model(
    file::String,
    model_type::Type,
    optimizer::Any,
    build_method::Function;
    kwargs...,
)::Dict{String,<:Any}
    return solve_model(parse_file(file), model_type, optimizer, build_method; kwargs...)
end


"""
    solve_model(
        data::Dict{String,<:Any},
        model_type::Type,
        optimizer::Any,
        build_method::Function;
        ref_extensions::Vector{<:Function}=Vector{Function}([]),
        solution_processors::Vector{<:Function}=Vector{Function}([]),
        relax_integrality::Bool=false,
        multinetwork::Bool=false,
        kwargs...)::Dict{String,<:Any}

Instantiates and solves the modeling object from data dictionary `data`. Here,
`model_type` is the model formulation type (e.g., NCWaterModel), `optimizer` is the
optimization solver used to solve the problem (e.g., Gurobi.Optimizer), and
`build_method` is the function used for building the problem specification being
considered (e.g., build_mn_owf). Moreover, `ref_extensions` is a vector of functions that
modify `ref`, `solution_processors` is a vector of functions that post-process model
solutions, `relax_integrality` is a Boolean indicating if the model solved should be
continuously relaxed, and `multinetwork` is a Boolean indicating whether or not the model
being solved is a multinetwork model. Returns a dictionary of optimization results.
"""
function solve_model(
    data::Dict{String,<:Any},
    model_type::Type,
    optimizer::Any,
    build_method::Function;
    ref_extensions::Vector{<:Function}=Vector{Function}([]),
    solution_processors::Vector{<:Function}=Vector{Function}([]),
    relax_integrality::Bool=false,
    multinetwork::Bool=false,
    kwargs...,
)::Dict{String,<:Any}
    if multinetwork != _IM.ismultinetwork(data)
        model_requirement = multinetwork ? "multi-network" : "single-network"
        data_type = _IM.ismultinetwork(data) ? "multi-network" : "single-network"
        message = "Attempted to build a $(model_requirement) model with $(data_type) data."
        Memento.error(_LOGGER, message)
    end

    wm = instantiate_model(
        data,
        model_type,
        build_method;
        ref_extensions=ref_extensions,
        kwargs...,
    )

    return optimize_model!(
        wm,
        optimizer=optimizer,
        relax_integrality=relax_integrality,
        solution_processors=solution_processors,
    )
end


"""
    instantiate_model(file::String, model_type::Type, build_method::Function; kwargs...)::AbstractWaterModel

Instantiates and returns a WaterModels object from data dictionary `data`. Here,
`model_type` is the formulation type (e.g., NCWaterModel) and `build_method` is the
function for building the problem specification being considered (e.g., build_mn_owf).
"""
function instantiate_model(
    file::String,
    model_type::Type,
    build_method::Function;
    kwargs...,
)::AbstractWaterModel
    return instantiate_model(parse_file(file), model_type, build_method; kwargs...)
end


"""
    instantiate_model(data::Dict{String,<:Any}, model_type::Type, build_method::Function; kwargs...)::AbstractWaterModel

Instantiates and returns a WaterModels object from input file with file path `file`.
Here, `model_type` is the formulation type (e.g., NCWaterModel) and `build_method` is the
function for building the problem specification being considered (e.g., build_mn_owf).
"""
function instantiate_model(
    data::Dict{String,<:Any},
    model_type::Type,
    build_method::Function;
    kwargs...,
)::AbstractWaterModel
    return _IM.instantiate_model(
        data, model_type, build_method, ref_add_core!, _wm_global_keys, wm_it_sym; kwargs...)
end


"""
Builds the ref dictionary from the data dictionary. Additionally, the ref dictionary can
contain the fields populated by the optional vector of `ref_extensions` functions.
"""
function build_ref(
    data::Dict{String,<:Any};
    ref_extensions::Vector{<:Function}=Vector{Function}([]),
)
    return _IM.build_ref(
        data,
        ref_add_core!,
        _wm_global_keys,
        wm_it_name,
        ref_extensions=ref_extensions,
    )
end


"""
    ref_add_core!(ref::Dict{Symbol,<:Any}; kwargs...)

Populates the dictionary `ref`, which stores commonly-used, precomputed data. The
function also converts data types, filters deactivated components, and stores values that
are computed from systemwide analyses of the original input data. Some of the most
important keys of this dictionary describe common network components, including:
* `:pipe` -- the set of pipes,
* `:des_pipe` -- the set of design pipes,
* `:pump` -- the set of pumps,
* `:regulator` -- the set of pressure regulating valves,
* `:short_pipe` -- the set of short pipes,
* `:valve` -- the set of gate and check valves,
* `:node` -- the set of nodes,
* `:demand` -- the set of demands,
* `:reservoir` -- the set of reservoirs,
* `:tank` -- the set of tanks
"""
function ref_add_core!(ref::Dict{Symbol,<:Any})
    _ref_add_core!(ref[:it][wm_it_sym][:nw], ref[:it][wm_it_sym][:head_loss])
end


function _filter_active_components(components::Dict{Int,<:Any})::Dict{Int,<:Any}
    return filter(x -> x.second["status"] !== STATUS_INACTIVE, components)
end


function _build_node_map(nodes::Dict{Int,<:Any}, components::Dict{Int,<:Any})
    ref_fr = Dict{Int,Array{Int,1}}(i => Array{Int,1}([]) for i in keys(nodes))
    ref_to = Dict{Int,Array{Int,1}}(i => Array{Int,1}([]) for i in keys(nodes))

    for (i, component) in components
        push!(ref_fr[component["node_fr"]], i)
        push!(ref_to[component["node_to"]], i)
    end

    return ref_fr, ref_to
end


function _pumps_match(pump_1::Dict{String, <:Any}, pump_2::Dict{String, <:Any})
    if sort(collect(keys(pump_1))) != sort(collect(keys(pump_2)))
        return false
    else
        # TODO: We should be more careful here.
        for key in ["node_fr", "node_to"]
            if pump_1[key] != pump_2[key]
                return false
            end
        end
    end

    return true
end


function _build_pump_groups(pumps::Dict{Int, <:Any})
    pump_group_indices = Set([])

    for (i, pump) in pumps
        other_pumps = filter(x -> x.first != i, pumps)
        matching_pumps = filter(x -> _pumps_match(pump, x.second), other_pumps)

        if length(matching_pumps) > 0
            pump_indices = sort(collect(vcat(i, keys(matching_pumps)...)))
            push!(pump_group_indices, Set(pump_indices))
        end
    end

    pump_group_indices = collect(pump_group_indices)
    return Dict{Int, Any}(i => Dict{String, Any}("pump_indices" =>
        pump_group_indices[i]) for i in 1:length(pump_group_indices))
end


function _set_ref_pump_head_gain_properties!(pumps::Dict{Int, <:Any})
    map(x -> x["head_curve_function"] = _calc_head_curve_function(x), values(pumps))
    map(x -> x["head_curve_derivative"] = _calc_head_curve_derivative(x), values(pumps))
    map(x -> x["head_curve_coefficients"] = _calc_head_curve_coefficients(x), values(pumps))
end


function _ref_add_core!(nw_refs::Dict{Int,<:Any}, head_loss::String)
    for (nw, ref) in nw_refs
        # Remove inactive nodes from the ref data dictionary.
        ref[:node] = _filter_active_components(ref[:node])

        # Create mappings of "from" and "to" arcs for link- (i.e., edge-) type components.
        for name in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
            ref[Symbol(name)] = _filter_active_components(ref[Symbol(name)])
            fr_sym, to_sym = Symbol(name * "_fr"), Symbol(name * "_to")
            ref[fr_sym], ref[to_sym] = _build_node_map(ref[:node], ref[Symbol(name)])
        end

        # Collect common arcs (i.e., node pairs) of design pipes in the network
        des_arcs = collect(Set((x["node_fr"], x["node_to"]) for (i, x) in ref[:des_pipe]))
        ref[:des_pipe_arc] = Dict{Int,Any}(i => des_arcs[i] for i in 1:length(des_arcs))
        
        # Set up dictionaries mapping node indices to attached component indices.
        for name in ["demand", "reservoir", "tank"]
            # Filter inactive components from the ref data dictionary.
            ref[Symbol(name)] = _filter_active_components(ref[Symbol(name)])

            # Initialize a dictionary to store a mapping of nodes to nodal components.
            name_sym = Symbol("node_" * name)
            ref[name_sym] = Dict{Int,Array{Int,1}}(i => Int[] for (i, node) in ref[:node])

            # Populate the dictionary mapping nodes to arrays of nodal components.
            for (i, comp) in ref[Symbol(name)]
                push!(ref[name_sym][comp["node"]], i)
            end
        end

        # Collect dispatchable and nondispatchable demands in the network.
        ref[:dispatchable_demand] = filter(x -> x.second["dispatchable"], ref[:demand])
        ref[:nondispatchable_demand] = filter(x -> !x.second["dispatchable"], ref[:demand])
        ref[:pump_group] = _build_pump_groups(ref[:pump])

        # Set pump head gain functions and derivatives.
        _set_ref_pump_head_gain_properties!(ref[:pump])

        # Store the exponent used within head loss relationships.
        ref[:alpha] = uppercase(head_loss) == "H-W" ? 1.852 : 2.0
    end
end


# Helper functions for multinetwork AbstractWaterModel objects.
ismultinetwork(wm::AbstractWaterModel) = _IM.ismultinetwork(wm, wm_it_sym)
nw_ids(wm::AbstractWaterModel) = _IM.nw_ids(wm, wm_it_sym)
nws(wm::AbstractWaterModel) = _IM.nws(wm, wm_it_sym)


# Helper functions for AbstractWaterModel component indices.
ids(wm::AbstractWaterModel, nw::Int, key::Symbol) = _IM.ids(wm, wm_it_sym, nw, key)
ids(wm::AbstractWaterModel, key::Symbol; nw::Int=nw_id_default) = _IM.ids(wm, wm_it_sym, key; nw=nw)


# Helper functions for AbstractWaterModel `ref` access.
ref(wm::AbstractWaterModel, nw::Int=nw_id_default) = _IM.ref(wm, wm_it_sym, nw)
ref(wm::AbstractWaterModel, nw::Int, key::Symbol) = _IM.ref(wm, wm_it_sym, nw, key)
ref(wm::AbstractWaterModel, nw::Int, key::Symbol, idx) = _IM.ref(wm, wm_it_sym, nw, key, idx)
ref(wm::AbstractWaterModel, nw::Int, key::Symbol, idx, param::String) = _IM.ref(wm, wm_it_sym, nw, key, idx, param)
ref(wm::AbstractWaterModel, key::Symbol; nw::Int=nw_id_default) = _IM.ref(wm, wm_it_sym, key; nw=nw)
ref(wm::AbstractWaterModel, key::Symbol, idx; nw::Int=nw_id_default) = _IM.ref(wm, wm_it_sym, key, idx; nw=nw)
ref(wm::AbstractWaterModel, key::Symbol, idx, param::String; nw::Int=nw_id_default) = _IM.ref(wm, wm_it_sym, key, idx, param; nw=nw)


# Helper functions for AbstractWaterModel `var` access.
var(wm::AbstractWaterModel, nw::Int=nw_id_default) = _IM.var(wm, wm_it_sym, nw)
var(wm::AbstractWaterModel, nw::Int, key::Symbol) = _IM.var(wm, wm_it_sym, nw, key)
var(wm::AbstractWaterModel, nw::Int, key::Symbol, idx) = _IM.var(wm, wm_it_sym, nw, key, idx)
var(wm::AbstractWaterModel, key::Symbol; nw::Int=nw_id_default) = _IM.var(wm, wm_it_sym, key; nw=nw)
var(wm::AbstractWaterModel, key::Symbol, idx; nw::Int=nw_id_default) = _IM.var(wm, wm_it_sym, key, idx; nw=nw)


# Helper functions for AbstractWaterModel `con` access.
con(wm::AbstractWaterModel, nw::Int=nw_id_default) = _IM.con(wm, wm_it_sym; nw=nw)
con(wm::AbstractWaterModel, nw::Int, key::Symbol) = _IM.con(wm, wm_it_sym, nw, key)
con(wm::AbstractWaterModel, nw::Int, key::Symbol, idx) = _IM.con(wm, wm_it_sym, nw, key, idx)
con(wm::AbstractWaterModel, key::Symbol; nw::Int=nw_id_default) = _IM.con(wm, wm_it_sym, key; nw=nw)
con(wm::AbstractWaterModel, key::Symbol, idx; nw::Int=nw_id_default) = _IM.con(wm, wm_it_sym, key, idx; nw=nw)


# Helper functions for AbstractWaterModel `sol` access.
sol(wm::AbstractWaterModel, nw::Int, args...) = _IM.sol(wm, wm_it_sym, nw, args...)
sol(wm::AbstractWaterModel, args...; nw::Int=nw_id_default) = _IM.sol(wm, wm_it_sym, args...; nw=nw)
