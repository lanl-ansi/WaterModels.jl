# Constraint templates help simplify data wrangling across multiple Water
# Flow formulations by providing an abstraction layer between the network data
# and network constraint definitions.  The constraint template's job is to
# extract the required parameters from a given network data structure and
# pass the data as named arguments to the Water Flow formulations.
#
# Constraint templates should always be defined over "GenericWaterModel"
# and should never refer to model variables

"Constraint to ensure at least one direction is set to take flow away from a source."
function constraint_source_flow(wm::GenericWaterModel{T}, i::Int, n::Int = wm.cnw) where T
    connections = wm.ref[:nw][n][:connection]
    f_branches = collect(keys(filter(is_out_node_function(i), connections)))
    t_branches = collect(keys(filter(is_in_node_function(i), connections)))
    constraint_source_flow(wm, i, f_branches, t_branches, n)
end

"Constraint to ensure at least one direction is set to take flow to a junction with demand."
function constraint_sink_flow(wm::GenericWaterModel{T}, i::Int, n::Int = wm.cnw) where T
    connections = wm.ref[:nw][n][:connection]
    f_branches = collect(keys(filter(is_out_node_function(i), connections)))
    t_branches = collect(keys(filter(is_in_node_function(i), connections)))
    constraint_sink_flow(wm, i, f_branches, t_branches, n)
end
