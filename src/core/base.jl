"Root of the water formulation type hierarchy."
abstract type AbstractWaterModel <: _IM.AbstractInfrastructureModel end

"A macro for adding the base WaterModels fields to a type definition."
_IM.@def wm_fields begin WaterModels.@im_fields end

""
function run_model(file::String, model_type::Type, optimizer, build_method; kwargs...)
    data = parse_file(file)
    return run_model(data, model_type, optimizer, build_method; kwargs...)
end

""
function run_model(
    data::Dict{String,<:Any}, model_type::Type, optimizer, build_method::Function; ref_extensions=[], solution_processors=[], multinetwork::Bool=false, kwargs...)
    if multinetwork != _IM.ismultinetwork(data)
        model_requirement = multinetwork ? "multi-network" : "single-network"
        data_type = _IM.ismultinetwork(data) ? "multi-network" : "single-network"
        Memento.error(_LOGGER, "Attempted to build a $(model_requirement) model with $(data_type) data.")
    end

    wm = instantiate_model(data, model_type, build_method; ref_extensions=ref_extensions, kwargs...)
    return optimize_model!(wm, optimizer=optimizer, solution_processors=solution_processors)
end

""
function instantiate_model(file::String, model_type::Type, build_method; kwargs...)
    data = parse_file(file)
    return instantiate_model(data, model_type, build_method; kwargs...)
end

""
function instantiate_model(data::Dict{String,<:Any}, model_type::Type, build_method; kwargs...)
    return _IM.instantiate_model(data, model_type, build_method, ref_add_core!, _wm_global_keys; kwargs...)
end

"""
Builds the ref dictionary from the data dictionary. Additionally the ref
dictionary would contain fields populated by the optional vector of
ref_extensions provided as a keyword argument.
"""
function build_ref(data::Dict{String,<:Any}; ref_extensions=[])
    return _IM.build_ref(data, ref_add_core!, _wm_global_keys, ref_extensions=ref_extensions)
end

"""
Returns a dict that stores commonly-used, precomputed data from of the data
dictionary, primarily for converting data types, filtering out deactivated
components, and storing system-wide values that need to be computed globally.

Some of the common keys include:

* `:pipe` -- the set of pipes in the network without valves,
* `:check_valve` -- the set of pipes in the network with check valves,
* `:shutoff_valve` -- the set of pipes in the network wiith shutoff valves,
* `:pressure_reducing_valve` -- the set of pressure reducing valves in the network,
* `:pump` -- the set of pumps in the network,
* `:node` -- the set of all nodes in the network,
* `:junction` -- the set of junctions in the network,
* `:reservoir` -- the set of reservoirs in the network,
* `:tank` -- the set of tanks in the network
"""
function ref_add_core!(refs::Dict{Symbol,<:Any})
    _ref_add_core!(refs[:nw])
end

function _ref_add_core!(nw_refs::Dict{Int,<:Any})
    for (nw, ref) in nw_refs
        ref[:resistance] = calc_resistances(ref[:pipe], ref[:viscosity], ref[:head_loss])
        ref[:resistance_cost] = calc_resistance_costs(ref[:pipe], ref[:viscosity], ref[:head_loss])

        # TODO: Check and shutoff valves should not have pipe properties.
        edge_index = length(ref[:pipe]) > 0 ? maximum(collect(keys(ref[:pipe]))) + 1 : 1
        ref[:check_valve] = filter(has_check_valve, ref[:pipe])
        ref[:shutoff_valve] = filter(has_shutoff_valve, ref[:pipe])
        ref[:pipe] = filter(!has_check_valve, ref[:pipe])
        ref[:pipe] = filter(!has_shutoff_valve, ref[:pipe])

        ref[:pressure_reducing_valve] = filter(is_pressure_reducing_valve, ref[:valve])
        ref[:pipe_des] = filter(is_des_pipe, ref[:pipe])
        ref[:pipe_fixed] = filter(!is_des_pipe, ref[:pipe])

        # Modify the network for standard modeling of tanks.
        for (i, tank) in ref[:tank]
            # Create a new node, which will be connected to the tank with a shutoff valve.
            node = deepcopy(ref[:node][tank["tank_node"]])
            node["index"] = sort(collect(keys(ref[:node])))[end] + 1
            node["source_id"] = ["node", "tank_$(i)_dummy"]
            ref[:node][node["index"]] = node

            # Instantiate the properties that define the shutoff valve.
            sv = Dict{String,Any}("flow_direction"=>0, "node_fr"=>tank["tank_node"],
                "node_to"=>node["index"], "diameter"=>1.0, "minor_loss"=>0.0, "length"=>0.0,
                "name"=>"tank_$(i)_sv", "status"=>"SV",
                "source_id"=>["shutoff_valve", "tank_$(i)_sv"], "initial_status"=>"Closed",
                "control"=>Dict{String,Any}(), "roughness"=>0.0)

            # Set the tank node index to the index of the dummy node.
            tank["tank_node"] = node["index"]

            # Give the new shutoff valve a new index.
            sv["index"] = edge_index

            # Update portions of `ref` with shutoff valve data.
            ref[:resistance][sv["index"]] = [0.0]
            ref[:resistance_cost][sv["index"]] = [0.0]
            ref[:shutoff_valve][sv["index"]] = sv

            edge_index += 1
        end

        # Create mappings of "from" and "to" arcs for link- (i.e., edge-) type components.
        for name in ["check_valve", "shutoff_valve", "pipe", "pressure_reducing_valve", "pump"]
            fr_sym, to_sym = Symbol(name * "_fr"), Symbol(name * "_to")
            ref[fr_sym] = [(i, c["node_fr"], c["node_to"]) for (i, c) in ref[Symbol(name)]]
            ref[to_sym] = [(i, c["node_to"], c["node_fr"]) for (i, c) in ref[Symbol(name)]]
        end

        # Set up dictionaries mapping node indices to attached component indices.
        for name in ["junction", "tank", "reservoir"]
            name_sym = Symbol("node_" * name)
            ref[name_sym] = Dict{Int,Array{Int,1}}(i=>Int[] for (i, node) in ref[:node])

            for (i, comp) in ref[Symbol(name)]
                push!(ref[name_sym][comp["$(name)_node"]], i)
            end
        end

        # Set alpha, that is, the exponent used for head loss relationships.
        ref[:alpha] = uppercase(ref[:head_loss]) == "H-W" ? 1.852 : 2.0
    end
end
