# Constraint templates help simplify data wrangling across multiple Water
# Flow formulations by providing an abstraction layer between the network data
# and network constraint definitions.  The constraint template's job is to
# extract the required parameters from a given network data structure and
# pass the data as named arguments to the Water Flow formulations.
#
# Constraint templates should always be defined over "GenericWaterModel"
# and should never refer to model variables

#function constraint_flow_conservation{T}(wm::GenericWaterModel{T}, i, n::Int = wm.cnw)
#    out_arcs = collect(keys(filter((id, pipe) -> pipe["node1"] == i, wm.ref[:nw][n][:pipes])))
#    out_dirs = [wm.data["pipes"][a]["flow_direction"] for a in out_arcs]
#    in_arcs = collect(keys(filter((id, pipe) -> pipe["node2"] == i, wm.ref[:nw][n][:pipes])))
#    in_dirs = [wm.data["pipes"][a]["flow_direction"] for a in in_arcs]
#
#    all_dirs = [out_dirs; in_dirs]
#    known_directions = all(y->y != UNKNOWN, all_dirs)
#
#    if known_directions
#        constraint_flow_conservation_known_directions(wm, i, n)
#    else
#        constraint_flow_conservation_unknown_directions(wm, i, n)
#    end
#end

#function constraint_potential_flow_coupling{T}(wm::GenericWaterModel{T}, a,
#                                               relaxed::Bool = true,
#                                               num_separators::Int = 5,
#                                               n::Int = wm.cnw)
#    headloss_type = wm.ref[:nw][n][:options]["headloss"]
#
#    # Apply the correct constraint based on the friction loss formulation.
#    if headloss_type == "h-w"
#        if relaxed
#            constraint_potential_flow_coupling_relaxed_hw(wm, a, num_separators, n)
#        else
#            constraint_potential_flow_coupling_exact_hw(wm, a, n)
#        end
#    elseif headloss_type == "d-w"
#        if relaxed
#            constraint_potential_flow_coupling_relaxed_dw(wm, a, num_separators, n)
#        else
#            constraint_potential_flow_coupling_exact_dw(wm, a, n)
#        end
#    end
#end
#
#function constraint_define_gamma{T}(wm::GenericWaterModel{T},
#                                    a, n::Int = wm.cnw)
#    flow_direction = wm.ref[:nw][n][:pipes][a]["flow_direction"]
#
#    # Apply the correct constraint based on the flow direction.
#    if flow_direction == UNKNOWN
#        constraint_define_gamma_unknown_direction(wm, a, n)
#    elseif flow_direction == POSITIVE
#        constraint_define_gamma_positive_direction(wm, a, n)
#    elseif flow_direction == NEGATIVE
#        constraint_define_gamma_negative_direction(wm, a, n)
#    end
#end
#
#function constraint_flow_direction{T}(wm::GenericWaterModel{T}, a,
#                                      n::Int = wm.cnw)
#    flow_direction = wm.ref[:nw][n][:pipes][a]["flow_direction"]
#
#    # Apply the correct constraint based on the flow direction.
#    if flow_direction == UNKNOWN
#        constraint_bidirectional_flow(wm, a, n)
#    elseif flow_direction == POSITIVE
#        constraint_positive_flow(wm, a, n)
#    elseif flow_direction == NEGATIVE
#        constraint_negative_flow(wm, a, n)
#    end
#end
