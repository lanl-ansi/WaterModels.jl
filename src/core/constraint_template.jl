# Constraint templates help simplify data wrangling across multiple Water Flow
# formulations by providing an abstraction layer between the network data and
# network constraint definitions. The constraint template's job is to extract
# the required parameters from a given network data structure and pass the data
# as named arguments to the Water Flow formulations.
#
# Constraint templates should always be defined over "GenericWaterModel" and
# should never refer to model variables.


### Node Constraints ###
function constraint_flow_conservation(wm::GenericWaterModel{T}, i::Int; n::Int=wm.cnw) where T <: AbstractDirectedFlowFormulation
    # Create the constraint dictionary if necessary.
    if !haskey(con(wm, n), :flow_conservation)
        con(wm, n)[:flow_conservation] = Dict{Int, JuMP.ConstraintRef}()
    end

    node_junctions = ref(wm, n, :node_junctions, i)
    node_arcs_fr = ref(wm, n, :node_arcs_fr, i)
    node_arcs_to = ref(wm, n, :node_arcs_to, i)
    node_reservoirs = ref(wm, n, :node_reservoirs, i)

    # TBD
    # for tid in ref(wm, n, :node_tanks, i)
    #     tank = ref(wm, n, :tanks, tid)
    #     # TODO add tank vars as loads
    # end

    node_demands = Dict(k => ref(wm, n, :junctions, k, "demand") for k in node_junctions)

    constraint_directed_flow_conservation(wm, n, i, node_arcs_fr, node_arcs_to, node_reservoirs, node_demands)
end


function constraint_flow_conservation(wm::GenericWaterModel{T}, i::Int; n::Int=wm.cnw) where T <: AbstractUndirectedFlowFormulation
    # Create the constraint dictionary if necessary.
    if !haskey(con(wm, n), :flow_conservation)
        con(wm, n)[:flow_conservation] = Dict{Int, JuMP.ConstraintRef}()
    end

    node_junctions = ref(wm, n, :node_junctions, i)
    node_arcs_fr = ref(wm, n, :node_arcs_fr, i)
    node_arcs_to = ref(wm, n, :node_arcs_to, i)
    node_reservoirs = ref(wm, n, :node_reservoirs, i)

    # TBD
    # for tid in ref(wm, n, :node_tanks, i)
    #     tank = ref(wm, n, :tanks, tid)
    #     # TODO add tank vars as loads
    # end

    node_demands = Dict(k => ref(wm, n, :junctions, k, "demand") for k in node_junctions)

    constraint_undirected_flow_conservation(wm, n, i, node_arcs_fr, node_arcs_to, node_reservoirs, node_demands)
end

#=
constraint_sink_flow(wm, i)
constraint_source_flow(wm, i)


### Junction Constraints ###


### Tank Constraints ###


### Reservoir Constraints ###



### Link Constraints ###
constraint_link_flow(wm, a)
constraint_link_flow_ne(wm, a)


### Pipe Constraints ###
constraint_potential_loss_pipe(wm, a)
constraint_potential_loss_pipe_ne(wm, a)
constraint_resistance_selection_ne(wm, a)

### Pump Constraints ###
constraint_potential_loss_pump(wm, a)
=#
