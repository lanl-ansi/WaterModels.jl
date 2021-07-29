WaterModels.jl Change Log
=========================

### v0.8.0
- Remove the convention of defining the _number_ of flow partitioning points for parameterizing linearized formulations.
- Let users define flow partitioning points for head loss and head gain functions using `flow_partition` entries within pipe, design pipe, and pump objects.
- Implement a `set_flow_partitions!` function that leverages PolyhedralRelaxations.jl to partition the space of flows for pipes, design pipes, and pumps.
- Implement per-unit conversions for network data and automatically compute per-unit transformations if not provided by the user.
- Add new global keys for per-unit conversion factors.
- Implement `pump_group` component sets that model symmetrical groups of pumps.
- Implement symmetry-breaking constraints for symmetrical groups of pumps.
- Refactor flow direction, pump model, and status enumerated types.
- Implement additional constraints to ensure tank level bounds will be satisfied after an Euler step, beginning from a fixed tank level.
- Implement constraints to model pump on and off switching limits.
- Implement new OWF specification that includes pump switching limits.
- Implement new functions for setting variable warm-start values.
- Implement `apply_wm!`, which wraps the InfrastructureModels `apply!` function.
- Remove scaled pump power variables and use direct power variables in per-unit, instead.
- Correctly distinguish between time _points_ and time _intervals_ used in multinetwork formulations.
- By default, do _not_ convert pipes that are short to short pipes.
- Refactor portions of the optimization-based bound tightening (OBBT) algorithm.
- Allow for usage of OBBT on multinetwork data sets.
- Add helper functions for relaxing select binary variables.

### v0.7.0
- Incorporate InfrastructureModels multi-infrastructure features.
- Rename `run_` methods to `solve_` and add deprecation warnings.
- Add pump models that obey linear power curves and quadratic head gain curves.
- Improve bounds for node-attached flow variables.

### v0.6.0
- Standardized network data conventions.
- Implemented `NCD` and `PWLRD` formulations.
- Removed `QRD` and `CQRD` formulations.
- Cleaned up remainder of the formulation hierarchy.
- Modified data conventions for components, especially design pipes.

### v0.5.0
- Rename formulation types.
- Decoupled `check_valve` and `shutoff_valve` objects from `pipe` objects.
- Introduced a new `valve` component.
- Introduced a new `short_pipe` component.
- Renamed `pressure_reducing_valve` to `regulator`.
- Decomposed component constraints a bit more.
- Implemented `QRD` and `CQRD` formulations.
- Refactored bound-tightening utility.
- Rename `junction` to `demand`.

### v0.4.0
- Removed unnecessary dependencies.
- Removed the `link` component abstraction.
- Removed complex control logic that was not being used.
- Refactored EPANET parser and corrected nodal attributes.
- Fixed various bugs related to component indices.
- Added new constraints to tighten `owf` formulations.
- Added utility function for "unbinarizing" indicator variables.
- Added utility function for optimization-based bound tightening.
- Simplified integration tests and removed previous tests.

### v0.3.0
- Updated to InfrastructureModels v0.5.
- Removed the Constrained Water Flow (`cwf`) problem specification.
- Renamed `post_` methods to `build_`.
- Renamed `ne` (network expansion) methods to `des` (design).
- Renamed `NCNC` formulation to `NC`.
- Completed LA implementations of constraints.
- Implemented pressure reducing valves.
- Implemented shutoff valves.
- Removed pump control (not operation) constraints from the Water flow (wf) problem specification.
- Migrating experimental formulations (CNC and MICP-E) to [WaterModelsAnnex.jl](https://github.com/lanl-ansi/WaterModelsAnnex.jl).

### v0.2.0
- Update to JuMP v0.21 and InfrastructureModels v0.4.
- Update to use JuMP/MOI status values.
- Change to a singular naming convention.
- Add upper bounds on package dependencies.
- Add draft Optimal Water Flow problem.
- Add draft LA formulations.
- Add draft LRD formulations.
- Add draft Water Flow (wf) problem formulations.
- Rename previous Water Flow (wf) problem formulations to Constrained Water Flow (cwf).

### v0.1.0
- Refactored nearly everything.
- Implemented initial CNC, MICP, and NCNC models.
- Began separating algorithmic components into [WaterModelsAnnex.jl](https://github.com/lanl-ansi/WaterModelsAnnex.jl).

### v0.0.1
- Implemented global optimization algorithm for network design problems constrained by Hazen-Williams head loss physics.
