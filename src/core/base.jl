# Functions commonly used in the construction of water models.

"Root of the water formulation type hierarchy."
abstract type AbstractWaterModel <: _IM.AbstractInfrastructureModel end

"A macro for adding the base WaterModels fields to a type definition."
_IM.@def wm_fields begin WaterModels.@im_fields end

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
function instantiate_model(data::Dict{String,<:Any}, model_type::Type, build_method; kwargs...)
    return _IM.instantiate_model(data, model_type, build_method, ref_add_core!, _wm_global_keys; kwargs...)
end

""
function build_ref(data::Dict{String,<:Any}; ref_extensions=[], kwargs...)
    return _IM.build_ref(data, ref_add_core!, _wm_global_keys, ref_extensions=ref_extensions)
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
function ref_add_core!(refs::Dict{Symbol,<:Any})
    _ref_add_core!(refs[:nw])
end

function _ref_add_core!(nw_refs::Dict{Int,<:Any})
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
