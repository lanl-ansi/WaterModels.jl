# Problem Specifications
## Water Flow (WF)

### Functions
```julia
# Create head loss functions, if necessary.
_function_head_loss(wm)
```

### Objective
```julia
objective_wf(wm)
```

### Variables
```julia
# Create head loss functions, if necessary.
variable_head(wm)
variable_flow(wm)
variable_pump_head_gain(wm)
variable_pump_power(wm)

# Indicator (status) variables.
variable_des_pipe_indicator(wm)
variable_pump_indicator(wm)
variable_regulator_indicator(wm)
variable_valve_indicator(wm)

# Create flow-related variables for node attachments.
variable_demand_flow(wm)
variable_reservoir_flow(wm)
variable_tank_flow(wm)
```

### Constraints
```julia
# Flow conservation at all nodes.
for (i, node) in ref(wm, :node)
    constraint_flow_conservation(wm, i)
    constraint_node_directionality(wm, i)
end

# Constraints on pipe flows, heads, and physics.
for (a, pipe) in ref(wm, :pipe)
    constraint_pipe_head(wm, a)
    constraint_pipe_head_loss(wm, a)
    constraint_pipe_flow(wm, a)
end

# Selection of design pipes along unique arcs.
for (k, arc) in ref(wm, :des_pipe_arc)
    constraint_des_pipe_flow(wm, k, arc[1], arc[2])
    constraint_des_pipe_head(wm, k, arc[1], arc[2])
    constraint_des_pipe_selection(wm, k, arc[1], arc[2])
end

# Constraints on design pipe flows, heads, and physics.
for (a, des_pipe) in ref(wm, :des_pipe)
    constraint_on_off_des_pipe_head(wm, a)
    constraint_on_off_des_pipe_head_loss(wm, a)
    constraint_on_off_des_pipe_flow(wm, a)
end

# Constraints on pump flows, heads, and physics.
for (a, pump) in ref(wm, :pump)
    constraint_on_off_pump_head(wm, a)
    constraint_on_off_pump_head_gain(wm, a)
    constraint_on_off_pump_flow(wm, a)
    constraint_on_off_pump_power(wm, a)
end

for (k, pump_group) in ref(wm, :pump_group)
    constraint_on_off_pump_group(wm, k)
end

# Constraints on short pipe flows and heads.
for (a, regulator) in ref(wm, :regulator)
    constraint_on_off_regulator_head(wm, a)
    constraint_on_off_regulator_flow(wm, a)
end

# Constraints on short pipe flows and heads.
for (a, short_pipe) in ref(wm, :short_pipe)
    constraint_short_pipe_head(wm, a)
    constraint_short_pipe_flow(wm, a)
end

# Constraints on tank volumes.
for (i, tank) in ref(wm, :tank)
    # Set the initial tank volume.
    constraint_tank_volume(wm, i)
end

# Constraints on valve flows and heads.
for (a, valve) in ref(wm, :valve)
    constraint_on_off_valve_head(wm, a)
    constraint_on_off_valve_flow(wm, a)
end
```

## Multinetwork Water Flow (MN WF)

### Functions
```julia
# Create head loss functions, if necessary.
_function_head_loss(wm)
```

### Objective
```julia
objective_wf(wm)
```

### Variables
```julia
# Get all network IDs in the multinetwork.
network_ids = sort(collect(nw_ids(wm)))
network_ids_inner = length(network_ids) > 1 ? network_ids[1:end-1] : network_ids

for n in network_ids_inner
    # Physical variables.
    variable_head(wm; nw=n)
    variable_flow(wm; nw=n)
    variable_pump_head_gain(wm; nw=n)
    variable_pump_power(wm; nw=n)

    # Indicator (status) variables.
    variable_des_pipe_indicator(wm; nw=n)
    variable_pump_indicator(wm; nw=n)
    variable_regulator_indicator(wm; nw=n)
    variable_valve_indicator(wm; nw=n)

    # Create flow-related variables for node attachments.
    variable_demand_flow(wm; nw=n)
    variable_reservoir_flow(wm; nw=n)
    variable_tank_flow(wm; nw=n)
end

if length(network_ids) > 1
    variable_head(wm; nw = network_ids[end])
end
```

### Constraints
```julia
# Get all network IDs in the multinetwork.
network_ids = sort(collect(nw_ids(wm)))
network_ids_inner = length(network_ids) > 1 ? network_ids[1:end-1] : network_ids

for n in network_ids_inner
    # Flow conservation at all nodes.
    for i in ids(wm, :node; nw=n)
        constraint_flow_conservation(wm, i; nw=n)
        constraint_node_directionality(wm, i; nw=n)
    end

    # Constraints on pipe flows, heads, and physics.
    for a in ids(wm, :pipe; nw=n)
        constraint_pipe_flow(wm, a; nw=n)
        constraint_pipe_head(wm, a; nw=n)
        constraint_pipe_head_loss(wm, a; nw=n)
    end

    # Constraints on pump flows, heads, and physics.
    for a in ids(wm, :pump; nw=n)
        constraint_on_off_pump_head(wm, a; nw=n)
        constraint_on_off_pump_head_gain(wm, a; nw=n)
        constraint_on_off_pump_flow(wm, a; nw=n)
        constraint_on_off_pump_power(wm, a; nw=n)
    end

    # Constraints on groups of parallel pumps.
    for k in ids(wm, :pump_group; nw=n)
        constraint_on_off_pump_group(wm, k; nw=n)
    end

    # Constraints on short pipe flows and heads.
    for a in ids(wm, :regulator; nw=n)
        constraint_on_off_regulator_head(wm, a; nw=n)
        constraint_on_off_regulator_flow(wm, a; nw=n)
    end

    # Constraints on short pipe flows and heads.
    for a in ids(wm, :short_pipe; nw=n)
        constraint_short_pipe_head(wm, a; nw=n)
        constraint_short_pipe_flow(wm, a; nw=n)
    end

    # Constraints on valve flows and heads.
    for a in ids(wm, :valve; nw=n)
        constraint_on_off_valve_head(wm, a; nw=n)
        constraint_on_off_valve_flow(wm, a; nw=n)
    end
end

# Start with the first network, representing the initial time step.
n_1 = network_ids[1]

# Constraints on tank volumes.
for i in ids(wm, :tank; nw = n_1)
    # Set initial conditions of tanks.
    constraint_tank_volume(wm, i; nw = n_1)
end

if length(network_ids) > 1
    # Constraints on tank volumes.
    for n_2 in network_ids[2:end]
        # Constrain tank volumes after the initial time index.
        for i in ids(wm, :tank; nw = n_2)
            constraint_tank_volume(wm, i, n_1, n_2)
        end

        # Update the first network used for integration.
        n_1 = n_2
    end
end
```

## Multinetwork Water Flow with Switching Constraints (MN WF SWITCHING)

### Inherited Functions, Variables, Constraints, and Objective
```julia
# Build the base multinetwork problem.
build_mn_wf(wm)
```

### Variables
```julia
# Get all network IDs in the multinetwork.
network_ids = sort(collect(nw_ids(wm)))

# Get the first network ID in the multinetwork.
n_1 = network_ids[1]

for n_2 in network_ids[2:end-1]
    # Add pump switching variables.
    variable_pump_switch_on(wm; nw = n_2)
    variable_pump_switch_off(wm; nw = n_2)

    n_1 = n_2
end
```

### Constraints
```julia
# Get all network IDs in the multinetwork.
network_ids = sort(collect(nw_ids(wm)))

# Get the first network ID in the multinetwork.
n_1 = network_ids[1]

for n_2 in network_ids[2:end-1]
    for a in ids(wm, :pump, nw = n_2)
        # Add constraints that define switching variables.
        constraint_pump_switch_on(wm, a, n_1, n_2)
        constraint_pump_switch_off(wm, a, n_1, n_2)
    end

    n_1 = n_2
end

for a in ids(wm, :pump; nw = network_ids[1])
    # Add constraints on the total number of pump switches.
    constraint_on_off_pump_switch(wm, a, network_ids[2:end-1])
end
```

## Optimal Water Flow (OWF)

### Inherited Functions, Variables, Constraints, and Objective
```julia
# Build the water flow problem.
build_wf(wm)
```

### Objective
```julia
# Add the optimal water flow objective.
objective_owf(wm)
```

## Multinetwork Optimal Water Flow (MN OWF)

### Inherited Functions, Variables, Constraints, and Objective
```julia
# Build the water flow problem.
build_mn_wf(wm)
```

### Objective
```julia
# Add the optimal water flow objective.
objective_owf(wm)
```

### Constraints
```julia
# Get all network IDs in the multinetwork.
network_ids = sort(collect(nw_ids(wm)))

# Ensure tanks recover their initial volume.
n_1, n_f = network_ids[1], network_ids[end]

for i in ids(wm, n_f, :tank)
    constraint_tank_volume_recovery(wm, i, n_1, n_f)
end
```

## Optimal Design (DES)

### Inherited Functions, Variables, Constraints, and Objective
```julia
# Build the water flow problem.
build_wf(wm)
```

### Objective
```julia
# Add the network design objective.
objective_des(wm)
```
