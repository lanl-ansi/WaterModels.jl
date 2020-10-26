WaterModels.jl Change Log
=========================

### v0.5.0
- Rename formulation types.
- Decoupled `check_valve` and `shutoff_valve` objects from `pipe` objects.
- Introduced a new `valve` component.
- Introduced a new `short_pipe` component.
- Renamed `pressure_reducing_valve` to `regulator`.
- Decomposed component constraints a bit more.
- Implemented `QRD` and `CQRD` formulations.
- Refactored bound-tightening utility.

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
