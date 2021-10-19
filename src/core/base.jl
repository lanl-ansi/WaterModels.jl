"Root of the water formulation type hierarchy."
abstract type AbstractWaterModel <: _IM.AbstractInfrastructureModel end


"Macro for adding base WaterModels fields to the type definition."
_IM.@def wm_fields begin
    WaterModels.@im_fields
end


"""
    solve_model(
        file::String,
        model_type::Type,
        optimizer::Any,
        build_method::Function;
        kwargs...)::Dict{String,<:Any}

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
    data = parse_file(file)
    return solve_model(data, model_type, optimizer, build_method; kwargs...)
end


"""
    solve_model(
        data::Dict{String,<:Any},
        model_type::Type,
        optimizer::Any,
        build_method::Function;
        ref_extensions::Vector{<:Function} = Vector{Function}([]),
        solution_processors::Vector{<:Function} = Vector{Function}([]),
        relax_integrality::Bool = false,
        multinetwork::Bool = false,
        kwargs...)::Dict{String, <:Any}

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
    ref_extensions::Vector{<:Function} = Vector{Function}([]),
    solution_processors::Vector{<:Function} = Vector{Function}([]),
    relax_integrality::Bool = false,
    multinetwork::Bool = false,
    kwargs...,
)::Dict{String,<:Any}
    # Ensure that, if a multinetwork is constructed, multinetwork data is being used.
    if multinetwork != _IM.ismultinetwork(data)
        model_requirement = multinetwork ? "multi-network" : "single-network"
        data_type = _IM.ismultinetwork(data) ? "multi-network" : "single-network"
        message = "Attempted to build a $(model_requirement) model with $(data_type) data."
        Memento.error(_LOGGER, message)
    end

    # Instantiate the WaterModels modeling object.
    @elapsed time_instantiate =
        wm = instantiate_model(
            data,
            model_type,
            build_method;
            ref_extensions = ref_extensions,
            kwargs...,
        )

    Memento.debug(_LOGGER, "Model build time: $(time_instantiate) seconds.")

    # Solve the model and postprocess its solution.
    @elapsed time_optimize =
        result = optimize_model!(
            wm,
            optimizer = optimizer,
            relax_integrality = relax_integrality,
            solution_processors = solution_processors,
        )

    Memento.debug(_LOGGER, "Model solution time: $(time_optimize) seconds.")

    # Return the optimization result.
    return result
end


"""
    instantiate_model(
        file::String,
        model_type::Type,
        build_method::Function;
        kwargs...)::AbstractWaterModel

Instantiates and returns a WaterModels object from data dictionary `data`. Here,
`model_type` is the formulation type (e.g., NCWaterModel) and `build_method` is the
function for building the problem specification being considered (e.g., build_owf).
"""
function instantiate_model(
    file::String,
    model_type::Type,
    build_method::Function;
    kwargs...,
)::AbstractWaterModel
    data = parse_file(file)
    return instantiate_model(data, model_type, build_method; kwargs...)
end


"""
    instantiate_model(
        data::Dict{String,<:Any},
        model_type::Type,
        build_method::Function;
        kwargs...)::AbstractWaterModel

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
        _wm_global_keys,
        wm_it_sym;
        kwargs...,
    )
end


"""
    build_ref(
        data::Dict{String,<:Any},
        ref_extensions::Vector{<:Function} = Vector{Function}([]))::Dict{Symbol,<:Any}

Builds the ref dictionary from the data dictionary. Additionally, the ref dictionary can
contain the fields populated by the optional vector of `ref_extensions` functions.
"""
function build_ref(
    data::Dict{String,<:Any};
    ref_extensions::Vector{<:Function} = Vector{Function}([]),
)::Dict{Symbol,<:Any}
    return _IM.build_ref(
        data,
        ref_add_core!,
        _wm_global_keys,
        wm_it_name,
        ref_extensions = ref_extensions,
    )
end


"""
    ref_add_core!(
        ref::Dict{Symbol,<:Any};
        kwargs...)

Populates the dictionary `ref`, which stores commonly-used, precomputed data. The
function also converts data types, filters deactivated components, and stores values that
are computed from systemwide analysis of the original input data. Some of the most
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


"Private version of `_ref_add_core!` that handles most of the `ref` building."
function _ref_add_core!(nw_refs::Dict{Int,<:Any}, head_loss::String)
    for ref in values(nw_refs)
        # Remove only active nodes from the ref data dictionary.
        ref[:node] = _filter_active_components(ref[:node])

        # Create mappings of "from" and "to" arcs for link- (i.e., edge-) type components.
        for name in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
            ref[Symbol(name)] = _filter_active_components(ref[Symbol(name)])
            fr_sym, to_sym = Symbol(name * "_fr"), Symbol(name * "_to")
            ref[fr_sym], ref[to_sym] = _build_node_map(ref[:node], ref[Symbol(name)])
        end

        # Collect common arcs (i.e., node pairs) of design pipes in the network
        des_arcs = collect(Set((x["node_fr"], x["node_to"]) for (i, x) in ref[:des_pipe]))
        ref[:des_pipe_arc] = Dict{Int,Any}(i => des_arcs[i] for i = 1:length(des_arcs))

        # Set up dictionaries mapping node indices to attached component indices.
        for name in ["demand", "reservoir", "tank"]
            # Filter only active components from the ref data dictionary.
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

        # Identify groups of pumps with identical properties.
        pump_groups_new = _build_pump_groups(ref[:pump])

        if haskey(ref, :pump_group)
            ref[:pump_group] = merge(ref[:pump_group], pump_groups_new)
        else
            ref[:pump_group] = pump_groups_new
        end

        # Set pump head gain functions and derivatives.
        _set_ref_pump_head_gain_properties!(ref[:pump])

        # Store the exponent used within head loss relationships.
        ref[:alpha] = _get_exponent_from_head_loss_form(head_loss)
    end
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


"Determine if two pumps have identical properties."
function _pumps_match(pump_1::Dict{String,<:Any}, pump_2::Dict{String,<:Any})::Bool
    if sort(collect(keys(pump_1))) != sort(collect(keys(pump_2)))
        return false
    else
        # TODO: More properties should be checked, here.
        for key in ["node_fr", "node_to"]
            if pump_1[key] != pump_2[key]
                return false
            end
        end
    end

    return true
end


"Compute groups of symmetric pumps with identical properties."
function _build_pump_groups(pumps::Dict{Int,<:Any})::Dict{Int,Any}
    # Initialize vectors of pump indices that can be grouped.
    pump_group_indices = Set([])

    for (i, pump) in pumps
        # Compute the pumps that have matching properties with `pump`.
        other_pumps = filter(x -> x.first != i, pumps)
        matching_pumps = filter(x -> _pumps_match(pump, x.second), other_pumps)

        if length(matching_pumps) > 0
            # If there are matching pumps, identify a pump group.
            pump_indices = sort(collect(vcat(i, keys(matching_pumps)...)))
            push!(pump_group_indices, Set(pump_indices))
        end
    end

    # Collect the vectors that describe all pump groups.
    pump_group_indices = collect([collect(x) for x in pump_group_indices])

    # Return the `:pump_group` object.
    return Dict{Int,Any}(
        i => Dict{String,Any}("pump_indices" => sort(pump_group_indices[i])) for
        i = 1:length(pump_group_indices)
    )
end


"Store pump head gain functional properties in each pump object."
function _set_ref_pump_head_gain_properties!(pumps::Dict{Int,<:Any})
    # Store head curve functional information within each pump object.
    map(x -> x["head_curve_function"] = _calc_head_curve_function(x), values(pumps))
    map(x -> x["head_curve_derivative"] = _calc_head_curve_derivative(x), values(pumps))
    map(x -> x["head_curve_coefficients"] = _calc_head_curve_coefficients(x), values(pumps))
end


# Helper functions for multinetwork AbstractWaterModel objects.
ismultinetwork(wm::AbstractWaterModel) = _IM.ismultinetwork(wm, wm_it_sym)
nw_ids(wm::AbstractWaterModel) = _IM.nw_ids(wm, wm_it_sym)
nws(wm::AbstractWaterModel) = _IM.nws(wm, wm_it_sym)


# Helper functions for AbstractWaterModel component indices.
ids(wm::AbstractWaterModel, nw::Int, key::Symbol) = _IM.ids(wm, wm_it_sym, nw, key)
ids(wm::AbstractWaterModel, key::Symbol; nw::Int = nw_id_default) =
    _IM.ids(wm, wm_it_sym, key; nw = nw)


# Helper functions for AbstractWaterModel `ref` access.
ref(wm::AbstractWaterModel, nw::Int = nw_id_default) = _IM.ref(wm, wm_it_sym, nw)
ref(wm::AbstractWaterModel, nw::Int, key::Symbol) = _IM.ref(wm, wm_it_sym, nw, key)
ref(wm::AbstractWaterModel, nw::Int, key::Symbol, idx) =
    _IM.ref(wm, wm_it_sym, nw, key, idx)
ref(wm::AbstractWaterModel, nw::Int, key::Symbol, idx, param::String) =
    _IM.ref(wm, wm_it_sym, nw, key, idx, param)
ref(wm::AbstractWaterModel, key::Symbol; nw::Int = nw_id_default) =
    _IM.ref(wm, wm_it_sym, key; nw = nw)
ref(wm::AbstractWaterModel, key::Symbol, idx; nw::Int = nw_id_default) =
    _IM.ref(wm, wm_it_sym, key, idx; nw = nw)
ref(wm::AbstractWaterModel, key::Symbol, idx, param::String; nw::Int = nw_id_default) =
    _IM.ref(wm, wm_it_sym, key, idx, param; nw = nw)


# Helper functions for AbstractWaterModel `var` access.
var(wm::AbstractWaterModel, nw::Int = nw_id_default) = _IM.var(wm, wm_it_sym, nw)
var(wm::AbstractWaterModel, nw::Int, key::Symbol) = _IM.var(wm, wm_it_sym, nw, key)
var(wm::AbstractWaterModel, nw::Int, key::Symbol, idx) =
    _IM.var(wm, wm_it_sym, nw, key, idx)
var(wm::AbstractWaterModel, key::Symbol; nw::Int = nw_id_default) =
    _IM.var(wm, wm_it_sym, key; nw = nw)
var(wm::AbstractWaterModel, key::Symbol, idx; nw::Int = nw_id_default) =
    _IM.var(wm, wm_it_sym, key, idx; nw = nw)


# Helper functions for AbstractWaterModel `con` access.
con(wm::AbstractWaterModel, nw::Int = nw_id_default) = _IM.con(wm, wm_it_sym; nw = nw)
con(wm::AbstractWaterModel, nw::Int, key::Symbol) = _IM.con(wm, wm_it_sym, nw, key)
con(wm::AbstractWaterModel, nw::Int, key::Symbol, idx) =
    _IM.con(wm, wm_it_sym, nw, key, idx)
con(wm::AbstractWaterModel, key::Symbol; nw::Int = nw_id_default) =
    _IM.con(wm, wm_it_sym, key; nw = nw)
con(wm::AbstractWaterModel, key::Symbol, idx; nw::Int = nw_id_default) =
    _IM.con(wm, wm_it_sym, key, idx; nw = nw)


# Helper functions for AbstractWaterModel `sol` access.
sol(wm::AbstractWaterModel, nw::Int, args...) = _IM.sol(wm, wm_it_sym, nw, args...)
sol(wm::AbstractWaterModel, args...; nw::Int = nw_id_default) =
    _IM.sol(wm, wm_it_sym, args...; nw = nw)
