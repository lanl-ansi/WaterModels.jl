# Constraints
```@meta
CurrentModule = WaterModels
```

## Constraint Templates
Constraint templates help simplify data wrangling across multiple optimization formulations by providing an abstraction layer between the network data and network constraint definitions.
Each constraint template's job is to extract the required parameters from a given network data structure and pass the data as arguments to the optimization formulation.
These templates should be defined over `AbstractWaterModel` and should not refer to model variables.
For more details, see the files `core/constraint_template.jl` and `core/constraint.jl`.
Here, `core/constraint_template.jl` provides higher-level constraint interfaces by processing network data and calling methods defined in `core/constraint.jl`.

## Nodal Constraints
```@docs
constraint_flow_conservation
constraint_sink_directionality
constraint_source_directionality
constraint_intermediate_directionality
```

## Tank Constraints
```@docs
constraint_tank_volume
constraint_tank_volume_fixed
constraint_tank_volume_recovery
```

## Pipe Constraints
```@docs
constraint_pipe_flow
constraint_pipe_head
constraint_pipe_head_loss
```

## Design Pipe Constraints
```@docs
constraint_des_pipe_flow
constraint_des_pipe_head
constraint_on_off_des_pipe_flow
constraint_on_off_des_pipe_head
constraint_on_off_des_pipe_head_loss
constraint_des_pipe_selection
```

## Short Pipe Constraints
```@docs
constraint_short_pipe_flow
constraint_short_pipe_head
```

## Pump Constraints
```@docs
constraint_on_off_pump_flow
constraint_on_off_pump_head
constraint_on_off_pump_head_gain
constraint_on_off_pump_power
constraint_on_off_pump_power_best_efficiency
constraint_on_off_pump_power_custom
constraint_on_off_pump_group
constraint_on_off_pump_switch
constraint_pump_switch_on
constraint_pump_switch_off
```

## Valve Constraints
```@docs
constraint_on_off_valve_flow
constraint_on_off_valve_head
```

## Regulator Constraints
```@docs
constraint_on_off_regulator_flow
constraint_on_off_regulator_head
```