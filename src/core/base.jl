"Root of the water formulation type hierarchy."
abstract type AbstractWaterModel end

"""
A macro for adding the base WaterModels fields to a type definition, where
* `data` is the original data, usually from reading in a `.inp` file,
* `setting` usually looks something like `Dict("output" => Dict("flows" => true))`, and
* `ref` is a place to store commonly-used precomputed data from the data dictionary,
  primarily for converting datatypes, filtering deactivated components, and storing
  system-wide values that need to be computed globally. See `build_ref(data)` for further details.
Methods on `AbstractWaterModel` for defining variables and adding constraints should
* work with the `ref` dict, rather than the original `data` dict,
* add them to `model::JuMP.AbstractModel`, and
* follow the conventions for variable and constraint names.
"""
InfrastructureModels.@def wm_fields begin
    model::JuMP.AbstractModel

    data::Dict{String,<:Any}
    setting::Dict{String,<:Any}
    solution::Dict{String,<:Any}

    ref::Dict{Symbol,<:Any}
    var::Dict{Symbol,<:Any}
    con::Dict{Symbol,<:Any}
    fun::Dict{Symbol,<:Any}
    cnw::Int

    # Extensions should define a type to hold information particular to
    # their functionality, and store an instance of the type in this
    # dictionary keyed on an extension-specific symbol.
    ext::Dict{Symbol,<:Any}
end

function InitializeWaterModel(WaterModel::Type, data::Dict{String,<:Any};
    ext=Dict{Symbol,Any}(), setting=Dict{String,Any}(),
    jump_model::JuMP.AbstractModel=JuMP.Model(), kwargs...)
    @assert WaterModel <: AbstractWaterModel

    ref = InfrastructureModels.ref_initialize(data, _wm_global_keys)

    var = Dict{Symbol, Any}(:nw => Dict{Int, Any}())
    con = Dict{Symbol, Any}(:nw => Dict{Int, Any}())
    fun = Dict{Symbol, Any}()

    for nw_id in keys(ref[:nw])
        var[:nw][nw_id] = Dict{Symbol, Any}()
        con[:nw][nw_id] = Dict{Symbol, Any}()
    end

    cnw = minimum([k for k in keys(var[:nw])])

    return WaterModel(
        jump_model,
        data,
        setting,
        Dict{String, Any}(), # Solution.
        ref,
        var,
        con,
        fun,
        cnw,
        ext
    )
end

# Helper functions for working with multinetworks.
""
ismultinetwork(wm::AbstractWaterModel) = (length(wm.ref[:nw]) > 1)

""
nw_ids(wm::AbstractWaterModel) = keys(wm.ref[:nw])

""
nws(wm::AbstractWaterModel) = wm.ref[:nw]

""
ids(wm::AbstractWaterModel, nw::Int, key::Symbol) = keys(wm.ref[:nw][nw][key])
ids(wm::AbstractWaterModel, key::Symbol; nw::Int=wm.cnw) = keys(wm.ref[:nw][nw][key])

""
ref(wm::AbstractWaterModel, nw::Int) = wm.ref[:nw][nw]
ref(wm::AbstractWaterModel, nw::Int, key::Symbol) = wm.ref[:nw][nw][key]
ref(wm::AbstractWaterModel, nw::Int, key::Symbol, idx) = wm.ref[:nw][nw][key][idx]
ref(wm::AbstractWaterModel, nw::Int, key::Symbol, idx, param::String) = wm.ref[:nw][nw][key][idx][param]

ref(wm::AbstractWaterModel; nw::Int=wm.cnw) = wm.ref[:nw][nw]
ref(wm::AbstractWaterModel, key::Symbol; nw::Int=wm.cnw) = wm.ref[:nw][nw][key]
ref(wm::AbstractWaterModel, key::Symbol, idx; nw::Int=wm.cnw) = wm.ref[:nw][nw][key][idx]


var(wm::AbstractWaterModel, nw::Int) = wm.var[:nw][nw]
var(wm::AbstractWaterModel, nw::Int, key::Symbol) = wm.var[:nw][nw][key]
var(wm::AbstractWaterModel, nw::Int, key::Symbol, idx) = wm.var[:nw][nw][key][idx]

var(wm::AbstractWaterModel; nw::Int=wm.cnw) = wm.var[:nw][nw]
var(wm::AbstractWaterModel, key::Symbol; nw::Int=wm.cnw) = wm.var[:nw][nw][key]
var(wm::AbstractWaterModel, key::Symbol, idx; nw::Int=wm.cnw) = wm.var[:nw][nw][key][idx]

""
con(wm::AbstractWaterModel, nw::Int) = wm.con[:nw][nw]
con(wm::AbstractWaterModel, nw::Int, key::Symbol) = wm.con[:nw][nw][key]
con(wm::AbstractWaterModel, nw::Int, key::Symbol, idx) = wm.con[:nw][nw][key][idx]

con(wm::AbstractWaterModel; nw::Int=wm.cnw) = wm.con[:nw][nw]
con(wm::AbstractWaterModel, key::Symbol; nw::Int=wm.cnw) = wm.con[:nw][nw][key]
con(wm::AbstractWaterModel, key::Symbol, idx; nw::Int=wm.cnw) = wm.con[:nw][nw][key][idx]

""
fun(wm::AbstractWaterModel) = wm.fun
fun(wm::AbstractWaterModel, key::Symbol) = wm.fun[key]

function JuMP.optimize!(wm::AbstractWaterModel, optimizer::JuMP.OptimizerFactory)
    if wm.model.moi_backend.state == _MOI.Utilities.NO_OPTIMIZER
        _, solve_time, solve_bytes_alloc, sec_in_gc = @timed JuMP.optimize!(wm.model, optimizer)
    else
        Memento.warn(_LOGGER, "Model already contains optimizer factory, cannot use optimizer specified in `solve_generic_model`")
        _, solve_time, solve_bytes_alloc, sec_in_gc = @timed JuMP.optimize!(wm.model)
    end

    try
        solve_time = _MOI.get(wm.model, _MOI.SolveTime())
    catch
        Memento.warn(_LOGGER, "the given optimizer does not provide the SolveTime() attribute, falling back on @timed.  This is not a rigorous timing value.");
    end

    return solve_time
end

""
function run_model(file::String, model_type::Type, optimizer, post_method; kwargs...)
    data = WaterModels.parse_file(file)
    return run_model(data, model_type, optimizer, post_method; kwargs...)
end

""
function run_model(data::Dict{String,<:Any}, model_type::Type, optimizer, post_method; ref_extensions=[], solution_builder=solution_owf!, kwargs...)
    #start_time = time()
    wm = build_model(data, model_type, post_method; ref_extensions=ref_extensions, kwargs...)
    #Memento.debug(_LOGGER, "wm model build time: $(time() - start_time)")

    #start_time = time()
    result = optimize_model!(wm, optimizer; solution_builder=solution_builder)
    #Memento.debug(_LOGGER, "wm model solve and solution time: $(time() - start_time)")

    return result
end

""
function build_model(file::String, model_type::Type, post_method; kwargs...)
    data = WaterModels.parse_file(file)
    return build_model(data, model_type, post_method; kwargs...)
end

""
function build_model(data::Dict{String,<:Any}, model_type::Type, post_method; ref_extensions=[], multinetwork=false, kwargs...)
    # NOTE, this model constructor will build the ref dict using the latest info from the data

    #start_time = time()
    wm = InitializeWaterModel(model_type, data; kwargs...)
    #Memento.info(LOGGER, "wm model_type time: $(time() - start_time)")

    if !multinetwork && ismultinetwork(wm)
        Memento.error(_LOGGER, "attempted to build a single-network model with multi-network data")
    end

    start_time = time()
    ref_add_core!(wm)
    for ref_ext in ref_extensions
        ref_ext(wm)
    end
    Memento.debug(_LOGGER, "wm build ref time: $(time() - start_time)")

    start_time = time()
    post_method(wm)
    Memento.debug(_LOGGER, "wm post_method time: $(time() - start_time)")

    return wm
end

""
function optimize_model!(wm::AbstractWaterModel, optimizer::JuMP.OptimizerFactory; solution_builder = solution_owf!)
    start_time = time()
    solve_time = JuMP.optimize!(wm, optimizer)
    Memento.debug(_LOGGER, "JuMP model optimize time: $(time() - start_time)")

    start_time = time()
    result = build_solution(wm, solve_time; solution_builder = solution_builder)
    Memento.debug(_LOGGER, "WaterModels solution build time: $(time() - start_time)")

    wm.solution = result["solution"]

    return result
end


"used for building ref without the need to build a initialize an AbstractWaterModel"
function build_ref(data::Dict{String,<:Any}; ref_extensions=[])
    ref = InfrastructureModels.ref_initialize(data, _wm_global_keys)
    _ref_add_core!(ref[:nw])

    for ref_ext in ref_extensions
        ref_ext(wm)
    end

    return ref
end

function ref_add_core!(wm::AbstractWaterModel)
    _ref_add_core!(wm.ref[:nw])
end

"""
Returns a dict that stores commonly-used, precomputed data from of the data
dictionary, primarily for converting data types, filtering out deactivated
components, and storing system-wide values that need to be computed globally.

Some of the common keys include:

* `:pipe` -- the set of pipes in the network,
* `:pump` -- the set of pumps in the network,
* `:valve` -- the set of valves in the network,
* `:link` -- the set of all links in the network,
* `:link_ne` -- the set of all network expansion links in the network,
* `:junction` -- the set of junctions in the network,
* `:reservoir` -- the set of reservoirs in the network,
* `:tank` -- the set of tanks in the network,
* `:emitters` -- the set of emitters in the network,
* `:node` -- the set of all nodes in the network
"""
function _ref_add_core!(nw_refs::Dict)
    for (nw, ref) in nw_refs
        ref[:link] = merge(ref[:pipe], ref[:valve], ref[:pump])
        ref[:pipe_ne] = filter(is_ne_link, ref[:pipe])
        ref[:link_ne] = filter(is_ne_link, ref[:link])
        ref[:check_valve] = filter(has_check_valve, ref[:pipe])

        # Set up arcs from existing links.
        ref[:arcs_fr] = [(i, comp["node_fr"], comp["node_to"]) for (i, comp) in ref[:link]]

        # Set up dictionaries mapping "from" links for a node.
        node_arcs_fr = Dict((i, Tuple{Int, Int, Int}[]) for (i, node) in ref[:node])

        for (l, i, j) in ref[:arcs_fr]
            push!(node_arcs_fr[i], (l, i, j))
        end

        ref[:node_arcs_fr] = node_arcs_fr

        # Set up dictionaries mapping "to" links for a node.
        node_arcs_to = Dict((i, Tuple{Int, Int, Int}[]) for (i, node) in ref[:node])

        for (l, i, j) in ref[:arcs_fr]
            push!(node_arcs_to[j], (l, i, j))
        end

        ref[:node_arcs_to] = node_arcs_to

        # Set up dictionaries mapping nodes to attached junctions.
        node_junctions = Dict((i, Int[]) for (i,node) in ref[:node])

        for (i, junction) in ref[:junction]
            push!(node_junctions[junction["junction_node"]], i)
        end

        ref[:node_junction] = node_junctions

        # Set up dictionaries mapping nodes to attached tanks.
        node_tanks = Dict((i, Int[]) for (i,node) in ref[:node])

        for (i,tank) in ref[:tank]
            push!(node_tanks[tank["tank_node"]], i)
        end

        ref[:node_tank] = node_tanks

        # Set up dictionaries mapping nodes to attached reservoirs.
        node_reservoirs = Dict((i, Int[]) for (i,node) in ref[:node])

        for (i,reservoir) in ref[:reservoir]
            push!(node_reservoirs[reservoir["reservoir_node"]], i)
        end

        ref[:node_reservoir] = node_reservoirs

        # TODO: Fix these when feeling ambitious about more carefully handling directions.
        #ref[:link_known_direction] = filter(has_known_flow_direction, ref[:link])
        #ref[:link_unknown_direction] = filter(!has_known_flow_direction, ref[:link])

        # Set the resistances based on the head loss type.
        headloss = ref[:option]["hydraulic"]["headloss"]
        viscosity = ref[:option]["hydraulic"]["viscosity"]
        ref[:resistance] = calc_resistances(ref[:pipe], viscosity, headloss)
        ref[:resistance_cost] = calc_resistance_costs(ref[:pipe], viscosity, headloss)

        # Set alpha, that is, the exponent used for head loss relationships.
        ref[:alpha] = uppercase(headloss) == "H-W" ? 1.852 : 2.0
    end
end

"Checks if any of the given keys are missing from the given dict."
function _check_missing_keys(dict, keys, formulation_type)
    missing_keys = []

    for key in keys
        if !haskey(dict, key)
            push!(missing_keys, key)
        end
    end

    if length(missing_keys) > 0
        error(_LOGGER, "The formulation $(formulation_type) requires the variable(s) $(keys), but the $(missing_keys) variable(s) were not found in the model.")
    end
end
