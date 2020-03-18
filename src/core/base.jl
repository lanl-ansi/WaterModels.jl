# Functions commonly used in the construction of all water models.

"Root of the water formulation type hierarchy."
abstract type AbstractWaterModel end

"""
A macro for adding the base WaterModels fields to a type definition, where
* `data` is the original data, usually from reading in a `.inp` file,
* `setting` usually looks something like `Dict("output"=>Dict("flows"=>true))`,
* `solution` is a dictionary to store model solutions with `data` conventions,
* `ref` is a dictionary to store commonly-used, precomputed data from the data
  dictionary, primarily for converting datatypes, filtering deactivated
  components, and storing system-wide values that need to be computed globally.
  See `build_ref(data)` for further details.
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

    sol::Dict{Symbol,<:Any}
    sol_proc::Dict{Symbol,<:Any}

    cnw::Int

    # Extensions should define a type to hold information particular to
    # their functionality, and store an instance of the type in this
    # dictionary keyed on an extension-specific symbol.
    ext::Dict{Symbol,<:Any}
end

"Default constructor for AbstractWaterModel types."
function InitializeWaterModel(WaterModel::Type, data::Dict{String,<:Any};
    ext=Dict{Symbol,Any}(), setting=Dict{String,Any}(),
    jump_model::JuMP.AbstractModel=JuMP.Model(), kwargs...)
    @assert WaterModel <: AbstractWaterModel

    ref = InfrastructureModels.ref_initialize(data, _wm_global_keys)
    var = Dict{Symbol, Any}(:nw=>Dict{Int,Any}())
    con = Dict{Symbol, Any}(:nw=>Dict{Int,Any}())
    fun = Dict{Symbol, Any}() # User-defined nonlinear functions.

    sol = Dict{Symbol,Any}(:nw=>Dict{Int,Any}())
    sol_proc = Dict{Symbol,Any}(:nw=>Dict{Int,Any}())

    for nw_id in keys(ref[:nw])
        nw_var = var[:nw][nw_id] = Dict{Symbol,Any}()
        nw_con = con[:nw][nw_id] = Dict{Symbol,Any}()
        nw_sol = sol[:nw][nw_id] = Dict{Symbol,Any}()
        nw_sol_proc = sol_proc[:nw][nw_id] = Dict{Symbol,Any}()
    end

    cnw = minimum([k for k in keys(var[:nw])])

    wm = WaterModel(
        jump_model,
        data,
        setting,
        Dict{String,Any}(), # solution
        ref,
        var,
        con,
        fun,
        sol,
        sol_proc,
        cnw,
        ext
    )

    return wm
end

report_duals(wm::AbstractWaterModel) = haskey(wm.setting, "output") &&
    haskey(wm.setting["output"], "duals") &&
    wm.setting["output"]["duals"] == true

# Helper functions for working with multinetworks.
"Returns whether or not the model is over a multinetwork."
ismultinetwork(wm::AbstractWaterModel) = (length(wm.ref[:nw]) > 1)

"Returns the network identifiers for multinetworks."
nw_ids(wm::AbstractWaterModel) = keys(wm.ref[:nw])

"Returns all multinetworks."
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
ref(wm::AbstractWaterModel, key::Symbol, idx, param::String; nw::Int=wm.cnw) = wm.ref[:nw][nw][key][idx][param]

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


""
sol(wm::AbstractWaterModel, nw::Int, args...) = _sol(wm.sol[:nw][nw], args...)
sol(wm::AbstractWaterModel, args...; nw::Int=wm.cnw) = _sol(wm.sol[:nw][nw], args...)

function _sol(sol::Dict, args...)
    for arg in args
        sol = haskey(sol, arg) ? sol[arg] : sol[arg] = Dict()
    end

    return sol
end

""
function run_model(file::String, model_type::Type, optimizer, build_method; kwargs...)
    data = WaterModels.parse_file(file)
    return run_model(data, model_type, optimizer, build_method; kwargs...)
end

""
function run_model(data::Dict{String,<:Any}, model_type::Type, optimizer,
    build_method; ref_extensions=[], solution_processors=[], kwargs...)
    wm = instantiate_model(data, model_type, build_method; ref_extensions=ref_extensions, kwargs...)
    return optimize_model!(wm, optimizer=optimizer, solution_processors=solution_processors)
end

""
function instantiate_model(file::String, model_type::Type, build_method; kwargs...)
    data = WaterModels.parse_file(file)
    return instantiate_model(data, model_type, build_method; kwargs...)
end

""
function instantiate_model(data::Dict{String,<:Any}, model_type::Type,
    build_method; ref_extensions=[], multinetwork=false, kwargs...)
    wm = InitializeWaterModel(model_type, data; kwargs...)

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
    build_method(wm)
    Memento.debug(_LOGGER, "wm build_method time: $(time() - start_time)")

    return wm
end

""
function optimize_model!(wm::AbstractWaterModel; optimizer=nothing, solution_processors=[])
    start_time = time()

    if optimizer != nothing
        if wm.model.moi_backend.state == _MOI.Utilities.NO_OPTIMIZER
            JuMP.set_optimizer(wm.model, optimizer)
        else
            Memento.warn(_LOGGER, "Model already contains optimizer, cannot "
                * "use optimizer specified in `optimize_model!`")
        end
    end

    if wm.model.moi_backend.state == _MOI.Utilities.NO_OPTIMIZER
        Memento.error(_LOGGER, "No optimizer specified in `optimize_model!` or the given JuMP model.")
    end

    _, solve_time, solve_bytes_alloc, sec_in_gc = @timed JuMP.optimize!(wm.model)

    try
        solve_time = _MOI.get(wm.model, _MOI.SolveTime())
    catch
        Memento.warn(_LOGGER, "The given optimizer does not provide the "
            * "SolveTime() attribute, falling back on @timed. This is not a "
            * "rigorous timing value.");
    end

    Memento.debug(_LOGGER, "JuMP model optimize time: $(time() - start_time)")
    start_time = time()
    result = build_result(wm, solve_time; solution_processors=solution_processors)
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
* `:emitter` -- the set of emitters in the network,
* `:node` -- the set of all nodes in the network
"""
function ref_add_core!(wm::AbstractWaterModel)
    _ref_add_core!(wm.ref[:nw])
end

function _ref_add_core!(nw_refs::Dict)
    for (nw, ref) in nw_refs
        ref[:link] = merge(ref[:pipe], ref[:valve], ref[:pump])
        ref[:check_valve] = filter(has_check_valve, ref[:pipe])
        ref[:link_ne] = filter(is_ne_link, ref[:link])
        ref[:pipe_ne] = filter(is_ne_link, ref[:pipe])
        ref[:link_fixed] = filter(!is_ne_link, ref[:link])
        ref[:pipe_fixed] = filter(!is_ne_link, ref[:pipe])

        # Set up arcs from existing links.
        ref[:arc_fr] = [(i, comp["node_fr"], comp["node_to"]) for (i, comp) in ref[:link]]

        # Set up dictionaries mapping "from" links for a node.
        node_arcs_fr = Dict((i, Tuple{Int, Int, Int}[]) for (i, node) in ref[:node])

        for (l, i, j) in ref[:arc_fr]
            push!(node_arcs_fr[i], (l, i, j))
        end

        ref[:node_arc_fr] = node_arcs_fr

        # Set up dictionaries mapping "to" links for a node.
        node_arcs_to = Dict((i, Tuple{Int, Int, Int}[]) for (i, node) in ref[:node])

        for (l, i, j) in ref[:arc_fr]
            push!(node_arcs_to[j], (l, i, j))
        end

        ref[:node_arc_to] = node_arcs_to

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
        error(_LOGGER, "The formulation $(formulation_type) requires the "
            * "variable(s) $(keys), but the $(missing_keys) variable(s) were "
            * "not found in the model.")
    end
end
