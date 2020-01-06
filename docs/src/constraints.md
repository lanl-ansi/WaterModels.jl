# Constraints
```@meta
CurrentModule = WaterModels
```

## Constraint Templates
Constraint templates help simplify data wrangling across multiple optimization formulations by providing an abstraction layer between the network data and network constraint definitions.
Each constraint template's job is to extract the required parameters from a given network data structure and pass the data as named arguments to the optimization formulation.

These templates should be defined over `AbstractWaterModel` and should not refer to model variables.
For more details, see the files `core/constraint_template.jl` and `core/constraint.jl`.
Here, `core/constraint_template.jl` provides higher-level constraint interfaces by processing network data and calling methods defined in `core/constraint.jl`.

## Nodal Constraints
```@docs
constraint_flow_conservation
```
