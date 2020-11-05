"Root of the water formulation type hierarchy."
abstract type AbstractWaterModel <: _IM.AbstractInfrastructureModel end


"A macro for adding the base WaterModels fields to a type definition."
_IM.@def wm_fields begin
    WaterModels.@im_fields
end


"""
    run_model(file, model_type, optimizer, build_method; kwargs...)

    Instantiates and solves the WaterModels modeling object from input file `file`. Here,
    `model_type` is the model formulation type, `optimizer` is the optimization solver used
    to solve the problem, and `build_method` is the build method for the problem
    specification being considered. Returns a dictionary of solution results.
"""
function run_model(
    file::String,
    model_type::Type,
    optimizer,
    build_method::Function;
    kwargs...,
)
    return run_model(parse_file(file), model_type, optimizer, build_method; kwargs...)
end


"""
    run_model(file, model_type, optimizer, build_method; kwargs...)

    Instantiates and solves the WaterModels modeling object from input data `data`. Here,
    `model_type` is the model formulation type, `optimizer` is the optimization solver used
    to solve the problem, and `build_method` is the build method for the problem
    specification being considered. Returns a dictionary of solution results.
"""
function run_model(
    data::Dict{String,<:Any},
    model_type::Type,
    optimizer,
    build_method::Function;
    ref_extensions::Vector{<:Function} = Vector{Function}([]),
    solution_processors::Vector{<:Function} = Vector{Function}([]),
    multinetwork::Bool = false,
    kwargs...,
)
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
        ref_extensions = ref_extensions,
        kwargs...,
    )

    return optimize_model!(
        wm,
        optimizer = optimizer,
        solution_processors = solution_processors,
    )
end


"""
    instantiate_model(file, model_type, build_method; kwargs...)

    Instantiates and returns WaterModels object from input file with path `file`. Here,
    `model_type` is the model formulation type and `build_method` is the build method for
    the problem specification being considered.
"""
function instantiate_model(
    file::String,
    model_type::Type,
    build_method::Function;
    kwargs...,
)
    return instantiate_model(parse_file(file), model_type, build_method; kwargs...)
end


"""
    instantiate_model(data, model_type, build_method; kwargs...)

    Instantiates and returns WaterModels object from input data `data`. Here, `model_type`
    is the model formulation type and `build_method` is the build method for the problem
    specification being considered.
"""
function instantiate_model(
    data::Dict{String,<:Any},
    model_type::Type,
    build_method::Function;
    kwargs...,
)
    return _IM.instantiate_model(
        data,
        model_type,
        build_method,
        ref_add_core!,
        _wm_global_keys;
        default_it = :wd,
        kwargs...,
    )
end


"""
Returns a dict that stores commonly-used, precomputed data from of the data dictionary,
primarily for converting data types, filtering out deactivated components, and storing
system-wide values that need to be computed globally.

Some of the common keys include:

* `:pipe` -- the set of pipes in the network,
* `:pump` -- the set of pumps in the network,
* `:regulator` -- the set of regulators in the network,
* `:short_pipe` -- the set of short pipes in the network,
* `:valve` -- the set of valves in the network,
* `:node` -- the set of all nodes in the network,
* `:demand` -- the set of demands in the network,
* `:reservoir` -- the set of reservoirs in the network,
* `:tank` -- the set of tanks in the network
"""
function ref_add_core!(refs::Dict{Symbol,<:Any})
    _ref_add_core!(refs[:it][:wd][:nw])
end


function _ref_add_core!(nw_refs::Dict{Int,<:Any})
    for (nw, ref) in nw_refs
        # Collect dispatchable and nondispatchable nodal components in the network.
        ref[:dispatchable_demand] = filter(x -> x.second["dispatchable"], ref[:demand])
        ref[:nondispatchable_demand] = filter(x -> !x.second["dispatchable"], ref[:demand])

        # Compute resistances for pipe-type components in the network.
        ref[:resistance] = calc_resistances(ref[:pipe], ref[:viscosity], ref[:head_loss])
        ref[:resistance_cost] =
            calc_resistance_costs(ref[:pipe], ref[:viscosity], ref[:head_loss])
        ref[:des_pipe] = filter(is_des_pipe, ref[:pipe])
        ref[:pipe] = filter(!is_des_pipe, ref[:pipe])

        # Create mappings of "from" and "to" arcs for link- (i.e., edge-) type components.
        for name in ["pipe", "pump", "regulator", "short_pipe", "valve"]
            fr_sym, to_sym = Symbol(name * "_fr"), Symbol(name * "_to")
            ref[fr_sym] = [(a, c["node_fr"], c["node_to"]) for (a, c) in ref[Symbol(name)]]
            ref[to_sym] = [(a, c["node_to"], c["node_fr"]) for (a, c) in ref[Symbol(name)]]
        end

        # Set up dictionaries mapping node indices to attached component indices.
        for name in ["demand", "tank", "reservoir"]
            name_sym = Symbol("node_" * name)
            ref[name_sym] = Dict{Int,Array{Int,1}}(i => Int[] for (i, node) in ref[:node])

            for (i, comp) in ref[Symbol(name)]
                push!(ref[name_sym][comp["node"]], i)
            end
        end

        # Set alpha, that is, the exponent used for head loss relationships.
        ref[:alpha] = uppercase(ref[:head_loss]) == "H-W" ? 1.852 : 2.0
    end
end
