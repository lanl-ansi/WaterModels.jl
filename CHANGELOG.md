WaterModels.jl Change Log
=========================

### v0.3.0
- Update to InfrastructureModels v0.4.3.
- Remove the Constrained Water Flow (cwf) problem specification.
- Rename `post_` methods to `solve_`.
- Rename `ne` (network expansion) methods to `des` (design).
- Complete MILP implementations of constraints.

### v0.2.0
- Update to JuMP v0.21 and InfrastructureModels v0.4.
- Update to use JuMP/MOI status values.
- Change to a singular naming convention.
- Add upper bounds on package dependencies.
- Add draft Optimal Water Flow problem.
- Add draft MILP formulations.
- Add draft MILP-R formulations.
- Add draft Water Flow (wf) problem formulations.
- Rename previous Water Flow (wf) problem formulations to Constrained Water Flow (cwf).

### v0.1.0
- Refactored nearly everything.
- Implemented initial CNLP, MICP, and NCNLP models.
- Began separating algorithmic components into [WaterModelsAnnex.jl](https://github.com/lanl-ansi/WaterModelsAnnex.jl).

### v0.0.1
- Implemented global optimization algorithm for network design problems constrained by Hazen-Williams head loss physics.
