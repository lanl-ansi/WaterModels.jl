# Functions for working with WaterModels data elements.


"""
Turns given single network data into multinetwork data with `count` replicates of the given
network. Note that this function performs a deepcopy of the network data. Significant
multinetwork space savings can often be achieved by building application specific methods of
building a multinetwork with minimal data replication (e.g., through storing references).
"""
function replicate(data::Dict{String,<:Any}, count::Int; global_keys::Set{String}=Set{String}())
    return _IM.replicate(data, count, union(global_keys, _wm_global_keys))
end


function _correct_flow_bounds!(data::Dict{String,<:Any})
    # TODO: Correct flow bounds for the remaining edge components, as well.
    for (idx, pump) in get(data, "pump", Dict{String,Any}())
        node_fr_id, node_to_id = string(pump["node_fr"]), string(pump["node_to"])
        node_fr, node_to = data["node"][node_fr_id], data["node"][node_to_id]
        pump["flow_min"] = _calc_pump_flow_min(pump, node_fr, node_to)
        pump["flow_max"] = _calc_pump_flow_max(pump, node_fr, node_to)
        pump["flow_min_forward"] = _calc_pump_flow_min_forward(pump, node_fr, node_to)
        pump["flow_max_reverse"] = _calc_pump_flow_max_reverse(pump, node_fr, node_to)
    end
end


"Turns a single network with a `time_series` data block into a multinetwork."
function make_multinetwork(data::Dict{String, <:Any}; global_keys::Set{String}=Set{String}())
    return InfrastructureModels.make_multinetwork(data, union(global_keys, _wm_global_keys))
end


function set_start!(data::Dict{String,<:Any}, component_type::String, var_name::String, start_name::String)
    if _IM.ismultinetwork(data)
        for (n, nw) in data["nw"]
            comps = values(nw[component_type])
            map(x -> x[start_name] = x[var_name], comps)
        end
    else
        comps = values(data[component_type])
        map(x -> x[start_name] = x[var_name], comps)
    end
end


function set_flow_start!(data::Dict{String,<:Any})
    set_start!(data, "pipe", "q", "q_pipe_start")
    set_start!(data, "pump", "q", "q_pump_start")
    set_start!(data, "regulator", "q", "q_regulator_start")
    set_start!(data, "short_pipe", "q", "q_short_pipe_start")
    set_start!(data, "valve", "q", "q_valve_start")
    set_start!(data, "reservoir", "q", "q_reservoir_start")
    set_start!(data, "tank", "q", "q_tank_start")
end


function set_head_start!(data::Dict{String,<:Any})
    set_start!(data, "node", "h", "h_start")
    set_start!(data, "pump", "g", "g_pump_start")
end


function set_indicator_start!(data::Dict{String,<:Any})
    set_start!(data, "pump", "status", "z_pump_start")
    set_start!(data, "regulator", "status", "z_regulator_start")
    set_start!(data, "valve", "status", "z_valve_start")
end


function set_start_all!(data::Dict{String,<:Any})
    set_flow_start!(data)
    set_head_start!(data)
    set_indicator_start!(data)
end


function fix_all_indicators!(data::Dict{String,<:Any})
    set_start!(data, "pump", "status", "z_min")
    set_start!(data, "pump", "status", "z_max")
    set_start!(data, "regulator", "status", "z_min")
    set_start!(data, "regulator", "status", "z_max")
    set_start!(data, "valve", "status", "z_min")
    set_start!(data, "valve", "status", "z_max")
end


"Translate a multinetwork dataset to a snapshot dataset with dispatchable components."
function _relax_network!(data::Dict{String,<:Any})
    _relax_nodes!(data)
    _relax_tanks!(data)
    _relax_reservoirs!(data)
    _relax_demands!(data)
end
