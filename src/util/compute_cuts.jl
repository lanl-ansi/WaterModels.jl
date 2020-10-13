function _get_binary_cut_mapping(data::Dict{String,<:Any}, optimizer; model_type::Type = CQRDWaterModel, ext::Dict{Symbol,<:Any} = Dict{Symbol,Any}())
    _make_reduced_data!(data)
    wm = instantiate_model(data, model_type, build_wf; ext=ext)
    indicator_index_set = _get_indicator_index_set(wm)
    JuMP.set_optimizer(wm.model, optimizer)
    binary_index_set = vcat(_get_indicator_index_set(wm), _get_direction_index_set(wm))
    binary_variable_mapping = []

    for id_1 in binary_index_set
        for id_2 in setdiff(binary_index_set, [id_1])
            match_0_min = _solve_coupled_binary(wm, id_1, id_2, 0.0, _MOI.MIN_SENSE)
            match_0_max = _solve_coupled_binary(wm, id_1, id_2, 0.0, _MOI.MAX_SENSE)
            match_1_min = _solve_coupled_binary(wm, id_1, id_2, 1.0, _MOI.MIN_SENSE)
            match_1_max = _solve_coupled_binary(wm, id_1, id_2, 1.0, _MOI.MAX_SENSE)

            if match_0_min == match_0_max
                append!(binary_variable_mapping, [((0.0, id_1), (match_0_min, id_2))])
            elseif match_1_min == match_1_max
                append!(binary_variable_mapping, [((1.0, id_1), (match_1_min, id_2))])
            end
        end
    end

    _revert_reduced_data!(data)
    return binary_variable_mapping
end


function _add_binary_cuts!(wm::AbstractWaterModel, binary_mapping)
    for nw in sort(collect(nw_ids(wm)))
        for mapping in binary_mapping
            index_1, index_2 = mapping[1][2], mapping[2][2]
            z_1 = var(wm, nw, index_1[3], index_1[4])
            z_2 = var(wm, nw, index_2[3], index_2[4])

            if mapping[1][1] == 0.0 && mapping[2][1] == 1.0
                JuMP.@constraint(wm.model, z_2 >= 1.0 - z_1)
            elseif mapping[1][1] == 1.0 && mapping[2][1] == 0.0
                JuMP.@constraint(wm.model, z_2 <= 1.0 - z_1)
            elseif mapping[1][1] == 1.0 && mapping[2][1] == 1.0
                JuMP.@constraint(wm.model, z_2 >= z_1)
            elseif mapping[1][1] == 0.0 && mapping[2][1] == 0.0
                JuMP.@constraint(wm.model, z_2 <= z_1)
            end
        end
    end
end


function _solve_coupled_binary(wm::AbstractWaterModel, index_1::Tuple, index_2::Tuple, fix_value::Float64, sense::_MOI.OptimizationSense)
    # Get the variable reference from the index tuple.
    z_1 = var(wm, index_1[1], index_1[3], index_1[4])
    z_2 = var(wm, index_2[1], index_2[3], index_2[4])

    # Fix the first variable.
    JuMP.fix(z_1, fix_value)

    # Optimize the variable (or affine expression) being tightened.
    JuMP.@objective(wm.model, sense, z_2)
    JuMP.optimize!(wm.model)
    termination_status = JuMP.termination_status(wm.model)

    if termination_status in [_MOI.LOCALLY_SOLVED, _MOI.OPTIMAL]
        candidate = JuMP.objective_value(wm.model) >= 0.5 ? 1.0 : 0.0
    else
        candidate = sense === _MOI.MIN_SENSE ? 0.0 : 1.0
    end

    # Unfix the first variable.
    JuMP.unfix(z_1)

    return candidate
end


function _sum_node_flows(wm::AbstractWaterModel, nw::Int64, i::Int64)
    # Collect various indices for node-type components connected to node `i`.
    reservoirs = ref(wm, nw, :node_reservoir, i) # Reservoirs attached to node `i`.
    tanks = ref(wm, nw, :node_tank, i) # Tanks attached to node `i`.
    demands = ref(wm, nw, :node_demand, i) # Demands attached to node `i`.

    # Sum the constant demands required at node `i`.
    nondispatchable_demands = filter(j -> j in ids(wm, nw, :nondispatchable_demand), demands)
    fixed_demands = [ref(wm, nw, :nondispatchable_demand, j)["flow_rate"] for j in nondispatchable_demands]
    net_fixed_demand = length(fixed_demands) > 0 ? sum(fixed_demands) : 0.0

    # Get the indices of dispatchable demands connected to node `i`.
    dispatchable_demands = filter(j -> j in ids(wm, nw, :dispatchable_demand), demands)

    # Generate the expression for flow.
    expr = JuMP.AffExpr(net_fixed_demand)
    expr -= length(tanks) > 0 ? sum(var(wm, nw, :q_tank, i) for i in tanks) : 0.0
    expr += length(dispatchable_demands) > 0 ? sum(var(wm, nw, :q_demand, i) for i in dispatchable_demands) : 0.0
    expr -= length(reservoirs) > 0 ? sum(var(wm, nw, :q_reservoir, i) for i in reservoirs) : 0.0

    return expr
end


function _add_flow_cut!(wm::AbstractWaterModel, n::Int64, var, nodes::Array{Int,1}, forward::Bool)
    flow_sum = JuMP.AffExpr(0.0)

    for i in nodes
        flow_sum += forward ? _sum_node_flows(wm, n, i) : -_sum_node_flows(wm, n, i)
    end

    JuMP.@constraint(wm.model, var == flow_sum)
end


function _add_flow_cuts!(wm::AbstractWaterModel)
    if _IM.ismultinetwork(wm.data)
        for (nw, nw_data) in wm.data["nw"]
            n = parse(Int64, nw)
            graph, node_map = WaterModels.create_graph(nw_data)
            reverse_node_map = Dict{Int,String}(i => x for (x, i) in node_map)

            for component_type in ["pipe", "pump", "regulator", "short_pipe", "valve"]
                for (a, comp) in nw_data[component_type]
                    i, j = node_map[string(comp["node_fr"])], node_map[string(comp["node_to"])]
                    LightGraphs.rem_edge!(graph, i, j)

                    if !LightGraphs.is_connected(graph)
                        graphs = LightGraphs.connected_components(graph)

                        for graph in filter(x -> length(x) > 1, graphs)
                            forward = i in graph ? false : true
                            q = var(wm, n, Symbol("q_" * component_type), comp["index"])
                            nodes = [parse(Int, reverse_node_map[i]) for i in graph]
                            _add_flow_cut!(wm, n, q, nodes, forward)
                        end
                    end

                    LightGraphs.add_edge!(graph, i, j)
                end
            end
        end
    else
        n = collect(nw_ids(wm))[1]
        graph, node_map = WaterModels.create_graph(wm.data)
        reverse_node_map = Dict{Int,String}(i => x for (x, i) in node_map)

        for component_type in ["pipe", "pump", "regulator", "short_pipe", "valve"]
            for (a, comp) in wm.data[component_type]
                i, j = node_map[string(comp["node_fr"])], node_map[string(comp["node_to"])]
                LightGraphs.rem_edge!(graph, i, j)

                if !LightGraphs.is_connected(graph)
                    graphs = LightGraphs.connected_components(graph)

                    for graph in filter(x -> length(x) > 1, graphs)
                        forward = i in graph ? false : true
                        q = var(wm, n, Symbol("q_" * component_type), comp["index"])
                        nodes = [parse(Int, reverse_node_map[i]) for i in graph]
                        _add_flow_cut!(wm, n, q, nodes, forward)
                    end
                end

                LightGraphs.add_edge!(graph, i, j)
            end
        end
    end
end
