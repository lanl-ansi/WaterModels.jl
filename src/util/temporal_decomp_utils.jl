########### functions ################
# function to build the hlinks dictionary
function build_temporal_hlink_dictionary(time_splits_array)
    path_length = size(time_splits_array)[1]
    hlinks_ids = 1:(path_length-1)
    hlinks_dict = Dict()
    for id in hlinks_ids
        hlinks_dict[id] = [id,id+1]
    end
    return hlinks_dict
end

function _build_node_subscriptions_dict(hlinks_dict)
    node_subscriptions_dict = Dict()
    for (hid,hlink) in hlinks_dict
        for nid in hlink
            if haskey(node_subscriptions_dict,nid) == true 
                push!(node_subscriptions_dict[nid],hid)
            else
                node_subscriptions_dict[nid] = []
                push!(node_subscriptions_dict[nid],hid)
            end
        end
    end
    return node_subscriptions_dict
end

function _build_node_consensus_dict(nid,hlink_subscriptions,hlinks_dict,wm_object)
    node_consenus_dict = Dict()

    array = collect(keys(wm_object.data["nw"]))
    times = [parse(Int, x) for x in array]
    for hid in hlink_subscriptions
        hlink = hlinks_dict[hid] 
        tuple_array = []
        if nid == hlink[1]
             consensus_idx = maximum(times)
        elseif nid == hlink[2]
            consensus_idx = minimum(times)
        else
            println("error")
        end
        tanks_dict = wm_object.data["nw"][string(consensus_idx)]["tank"]
        for (tank_id,tank) in tanks_dict
            node = tank["node"]
            pub_name = "$(consensus_idx)_h[$(node)]"
            private_name = "$(consensus_idx)_h[$(node)]"
            t = (pub_name,private_name)
            push!(tuple_array,t)
        end
        node_consenus_dict[hid] = tuple_array
    end

    return node_consenus_dict 
end


function build_temporal_consensus_dict(hlinks_dict,wms)
    node_subscriptions_dict = _build_node_subscriptions_dict(hlinks_dict)
    temporal_consensus_dict = Dict()
    for (nid,hlink_subscriptions) in node_subscriptions_dict
        wm_object = wms[nid]
        node_consenus_dict = _build_node_consensus_dict(nid,hlink_subscriptions,hlinks_dict,wm_object)
        temporal_consensus_dict[nid] = node_consenus_dict
    end
   return temporal_consensus_dict
end

