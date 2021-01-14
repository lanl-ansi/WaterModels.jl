"Root of the water formulation type hierarchy."
abstract type AbstractWaterModel <: _IM.AbstractInfrastructureModel end


"Macro for adding base WaterModels fields to the type definition."
_IM.@def wm_fields begin
    WaterModels.@im_fields
end


"""
    run_model(file::String, model_type::Type, optimizer::Any, build_method::Function; kwargs...)::Dict{String,<:Any}

Instantiates and solves the modeling object from input file with file path `file`. Here,
`model_type` is the model formulation type (e.g., NCWaterModel), `optimizer` is the
optimization solver used to solve the problem (e.g., Gurobi.Optimizer), and
`build_method` is the function used for building the problem specification being
considered (e.g., build_mn_owf). Returns a dictionary of optimization results.
"""
function run_model(
    file::String,
    model_type::Type,
    optimizer::Any,
    build_method::Function;
    kwargs...,
)::Dict{String,<:Any}
    return run_model(parse_file(file), model_type, optimizer, build_method; kwargs...)
end


"""
    run_model(data::Dict{String,<:Any}, model_type::Type, optimizer::Any,
              build_method::Function; ref_extensions::Vector{<:Function}=Vector{Function}([]),
              solution_processors::Vector{<:Function}=Vector{Function}([]),
              multinetwork::Bool=false, kwargs...)::Dict{String,<:Any}

Instantiates and solves the modeling object from data dictionary `data`. Here,
`model_type` is the model formulation type (e.g., NCWaterModel), `optimizer` is the
optimization solver used to solve the problem (e.g., Gurobi.Optimizer), and
`build_method` is the function used for building the problem specification being
considered (e.g., build_mn_owf). Moreover, `ref_extensions` is a vector of functions that
modify `ref`, `solution_processors` is a vector of functions that post-process model
solutions, and `multinetwork` is a Boolean indicating whether or not the model being
solved is a multinetwork model. Returns a dictionary of optimization results.
"""
function run_model(
    data::Dict{String,<:Any},
    model_type::Type,
    optimizer::Any,
    build_method::Function;
    ref_extensions::Vector{<:Function}=Vector{Function}([]),
    solution_processors::Vector{<:Function}=Vector{Function}([]),
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
        data,
        model_type,
        build_method,
        ref_add_core!,
        _wm_global_keys;
        kwargs...,
    )
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
    _ref_add_core!(ref[:nw], ref[:head_loss])
end


function _filter_active_components(components::Dict{Int,<:Any})::Dict{Int,<:Any}
    return filter(x -> x.second["status"] != 0, components)
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

        # Store the exponent used within head loss relationships.
        ref[:alpha] = uppercase(head_loss) == "H-W" ? 1.852 : 2.0
    end
end
