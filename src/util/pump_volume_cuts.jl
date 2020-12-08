function _collect_remaining_nws(wm::AbstractWaterModel, n::Int)
    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(nw_ids(wm)))
    n_id = findfirst(x -> x == n, network_ids)
    return network_ids[n_id:end]
end


function _sum_remaining_pump_flows(wm::AbstractWaterModel, nws::Array{Int64,1})
    return sum(sum(var(wm, nw, :q_pump)) for nw in nws)
end


function _sum_remaining_demands(wm::AbstractWaterModel, nws::Array{Int64,1})
    expr = JuMP.AffExpr(0.0)

    for nw in nws
        for (i, node) in ref(wm, nw, :node)
            # Sum the constant demands required at node `i`.
            demands = ref(wm, nw, :node_demand, i) # Demands attached to node `i`.
            nondispatchable_demands = filter(j -> j in ids(wm, nw, :nondispatchable_demand), demands)
            fixed_demands = [ref(wm, nw, :nondispatchable_demand, j)["flow_nominal"] for j in nondispatchable_demands]
            net_fixed_demand = length(fixed_demands) > 0 ? sum(fixed_demands) : 0.0

            # Get the indices of dispatchable demands connected to node `i`.
            dispatchable_demands = filter(j -> j in ids(wm, nw, :dispatchable_demand), demands)

            # Generate the expression for flow.
            expr += length(dispatchable_demands) > 0 ? sum(var(wm, nw, :q_demand, i) for i in dispatchable_demands) : 0.0
            expr += length(fixed_demands) > 0 ? sum(fixed_demands) : 0.0
        end
    end

    return expr
end


function _add_pump_volume_cuts!(wm::AbstractWaterModel)
    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(nw_ids(wm)))

    # Start with the first network, representing the initial time step.
    n_1 = network_ids[1]

    for (n, network) in nws(wm)
        nws_remaining = _collect_remaining_nws(wm, n)
        time_step = ref(wm, n, :time_step)
        demand_volume_remaining = _sum_remaining_demands(wm, nws_remaining) * time_step
        pump_volume_remaining = _sum_remaining_pump_flows(wm, nws_remaining) * time_step
        tank_volume_sum = sum(var(wm, n_1, :V)) - sum(var(wm, n, :V))
        c = JuMP.@constraint(wm.model, tank_volume_sum + demand_volume_remaining <= pump_volume_remaining)
    end
end
