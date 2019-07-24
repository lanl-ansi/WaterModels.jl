# Constraint templates help simplify data wrangling across multiple Water Flow
# formulations by providing an abstraction layer between the network data and
# network constraint definitions. The constraint template's job is to extract
# the required parameters from a given network data structure and pass the data
# as named arguments to the Water Flow formulations.
#
# Constraint templates should always be defined over "GenericWaterModel" and
# should never refer to model variables.


### Node Constraints ###
function constraint_flow_conservation(wm::GenericWaterModel, i::Int; nw::Int=wm.cnw)
    # Create the constraint dictionary if necessary.
    if !haskey(con(wm, nw), :flow_conservation)
        con(wm, nw)[:flow_conservation] = Dict{Int, JuMP.ConstraintRef}()
    end

    node_junctions = ref(wm, nw, :node_junctions, i)
    node_arcs_fr = ref(wm, nw, :node_arcs_fr, i)
    node_arcs_to = ref(wm, nw, :node_arcs_to, i)
    node_reservoirs = ref(wm, nw, :node_reservoirs, i)

    # TBD
    # for tid in ref(wm, n, :node_tanks, i)
    #     tank = ref(wm, n, :tanks, tid)
    #     # TODO add tank vars as loads
    # end

    node_demands = Dict(k => ref(wm, nw, :junctions, k, "demand") for k in node_junctions)

    constraint_flow_conservation(wm, nw, i, node_arcs_fr, node_arcs_to, node_reservoirs, node_demands)
end

#=
constraint_sink_flow(wm, i)
constraint_source_flow(wm, i)
=#

### Junction Constraints ###


### Tank Constraints ###


### Reservoir Constraints ###



### Link Constraints ###
function constraint_link_flow(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    if !haskey(con(wm, nw), :link_directed_flow)
        con(wm, nw)[:link_directed_flow] = Dict{Int, JuMP.ConstraintRef}()
    end

    constraint_link_flow(wm, nw, a)
end

function constraint_link_flow_ne(wm::GenericWaterModel, a::Int; nw::Int=wm.cnw)
    constraint_link_flow_ne(wm, nw, a)
end


#=
### Pipe Constraints ###
constraint_potential_loss_pipe(wm, a)
constraint_potential_loss_pipe_ne(wm, a)
constraint_resistance_selection_ne(wm, a)

### Pump Constraints ###
constraint_potential_loss_pump(wm, a)
=#
