# Problem Specifications
## Water Flow (WF)

### Functions
```julia
function_head_loss(wm)
```

### Objective
```julia
objective_wf(pm)
```

### Variables
```julia
variable_head(wm, bounded=false)
variable_flow(wm, bounded=false)
variable_volume(wm, bounded=false)
variable_reservoir(wm)
variable_tank(wm)
```

### Constraints
```julia
for (a, pipe) in ref(wm, :pipe)
    constraint_head_loss_pipe(wm, a)
end

for (a, pump) in ref(wm, :pump)
    constraint_head_gain_pump(wm, a, force_on=true)
end

for (i, node) in ref(wm, :node)
    constraint_flow_conservation(wm, i)
end

for (i, reservoir) in ref(wm, :reservoir)
    constraint_source_head(wm, i)
    constraint_source_flow(wm, i)
end

for (i, junction) in ref(wm, :junction)
    if junction["demand"] > 0.0
        constraint_sink_flow(wm, i)
    end
end

for i in ids(wm, :tank)
    constraint_link_volume(wm, i)
    constraint_tank_state(wm, i)
end
```
