# WaterModels.jl

## A Julia Package for Water Network Optimization

## Build Status
| [Linux][ci-link]  | [macOS][ci-link]  | [Codecov][cov-link]   |
| :---------------: | :---------------: | :-------------------: |
| ![ci-badge]       | ![ci-badge]       | ![cov-badge]          |

[ci-badge]: https://travis-ci.org/lanl-ansi/WaterModels.jl.svg?branch=master "Travis build status"
[ci-link]: https://travis-ci.org/lanl-ansi/WaterModels.jl "Travis build status"
[cov-badge]: https://codecov.io/gh/lanl-ansi/WaterModels.jl/branch/master/graph/badge.svg
[cov-link]: https://codecov.io/gh/lanl-ansi/WaterModels.jl

## Introduction
WaterModels.jl is a Julia package for steady state water network optimization.
It is designed to enable computational evaluation of historical and emerging water network formulations and algorithms using a common platform.
The software is engineered to decouple problem specifications (e.g., water flow, network expansion) from water network optimization formulations (e.g., mixed-integer linear, mixed-integer nonlinear).
This decoupling enables the definition of a wide variety of water network optimization formulations and their comparison on common problem specifications.

**Core Problem Specifications**
* Water Flow (wf)
* Network Expansion (ne)

**Core Network Formulations**
* CNLP (convex nonlinear program used for determining network flow rates)
* MICP (relaxation-based mixed-integer convex program)
* MILP-R (relaxation-based mixed-integer linear program)

## Usage at a Glance
To solve a network design problem constrained by the Hazen-Williams head loss relationship to global optimality, execute the following (using the `shamir` network as an example):
```
using GLPK
using GLPKMathProgInterface
using Ipopt
using WaterModels

glpk = GLPKSolverMIP(presolve = false, msg_lev = GLPK.MSG_OFF)
ipopt = IpoptSolver(print_level = 0, tol = 1.0e-9, max_iter=9999)

network_path = "test/data/epanet/shamir.inp"
modification_path = "test/data/json/shamir.json"
status = solve_global(network_path, modification_path, ipopt, glpk)
```
Using [Gurobi](https://github.com/JuliaOpt/Gurobi.jl) in place of [GLPK](https://github.com/JuliaOpt/GLPK.jl) should result in substantially faster convergence.

## Development
Community-driven development and enhancement of WaterModels is welcomed and encouraged.
Please feel free to fork this repository and share your contributions to the master branch with a pull request.
That said, it is important to keep in mind the code quality requirements and scope of WaterModels before preparing a contribution.
See [CONTRIBUTING.md](https://github.com/lanl-ansi/WaterModels.jl/blob/master/CONTRIBUTING.md) for code contribution guidelines.

## Acknowledgments
This code has been developed as part of the Advanced Network Science Initiative at Los Alamos National Laboratory.
The primary developer is [Byron Tasseff](https://github.com/tasseff) with support from the following contributors:
- [Russell Bent](https://github.com/rb004f), Los Alamos National Laboratory
- [Carleton Coffrin](https://github.com/ccoffrin), Los Alamos National Laboratory

## License
This code is provided under a [modified BSD license](https://github.com/lanl-ansi/WaterModels.jl/blob/master/LICENSE.md) as part of the Multi-Infrastructure Control and Optimization Toolkit (MICOT), LA-CC-13-108.
