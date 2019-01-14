# WaterModels.jl

## A Julia Package for Water Network Optimization

## Build Status
| [Linux][ci-link]  | [OSX][ci-link]    | [Codecov][cov-link]   |
| :---------------: | :---------------: | :-------------------: |
| ![ci-badge]       | ![ci-badge]       | ![cov-badge]          |

[ci-badge]: https://travis-ci.org/lanl-ansi/WaterModels.jl.svg?branch=master "Travis build status"
[ci-link]: https://travis-ci.org/lanl-ansi/WaterModels.jl "Travis build status"
[cov-badge]: https://codecov.io/gh/lanl-ansi/WaterModels.jl/branch/master/graph/badge.svg
[cov-link]: https://codecov.io/gh/lanl-ansi/WaterModels.jl

## Introduction
WaterModels.jl is a Julia package for steady state water network optimization.
It is designed to enable computational evaluation of emerging water network formulations and algorithms in a common platform.
The code is engineered to decouple problem specifications (e.g., feasibility, network expansion) from water network formulations (e.g., mixed-integer linear, mixed-integer nonlinear).
This enables the definition of a wide variety of water network optimization formulations and their comparison on common problem specifications.

**Core Problem Specifications**
* Water Flow (wf)
* Network Expansion (ne)

**Core Network Formulations**
* MINLP-R (relaxation-based mixed-integer nonlinear program)
* MILP-R (relaxation-based mixed-integer linear program)

## Highlights in v0.0.1 (2018-09-26)

## Usage at a Glance

## Development
Community-driven development and enhancement of WaterModels is welcomed and encouraged.
Please feel free to fork this repository and share your contributions to the master branch with a pull request.
That said, it is important to keep in mind the code quality requirements and scope of WaterModels before preparing a contribution.
See [CONTRIBUTING.md](https://github.com/lanl-ansi/WaterModels.jl/blob/master/CONTRIBUTING.md) for code contribution guidelines.

## Acknowledgments
This code has been developed as part of the Advanced Network Science Initiative at Los Alamos National Laboratory.
The primary developer is [Byron Tasseff](@tasseff) with support from the following contributors:
- [Russell Bent](@rb004f), Los Alamos National Laboratory
- [Carleton Coffrin](@ccoffrin), Los Alamos National Laboratory

## License
This code is provided under a [modified BSD license](https://github.com/lanl-ansi/WaterModels.jl/blob/master/LICENSE.md) as part of the Multi-Infrastructure Control and Optimization Toolkit (MICOT), LA-CC-13-108.
