# WaterModels.jl

<img src="https://lanl-ansi.github.io/WaterModels.jl/dev/assets/logo.svg" align="left" width="200" alt="WaterModels logo">

Release: [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://lanl-ansi.github.io/WaterModels.jl/stable)

Dev:
[![Build Status](https://travis-ci.org/lanl-ansi/WaterModels.jl.svg?branch=master)](https://travis-ci.org/lanl-ansi/WaterModels.jl)
[![codecov](https://codecov.io/gh/lanl-ansi/WaterModels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/lanl-ansi/WaterModels.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://lanl-ansi.github.io/WaterModels.jl/latest)

WaterModels.jl is a Julia package for steady state water network optimization.
It is designed to enable computational evaluation of historical and emerging water network formulations and algorithms using a common platform.
The software is engineered to decouple problem specifications (e.g., water flow, network expansion) from water network optimization formulations (e.g., mixed-integer linear, mixed-integer nonlinear).
This decoupling enables the definition of a wide variety of water network optimization formulations and their comparison on common problem specifications.

**Core Problem Specifications**
* Water Flow (wf)
* Constrained Water Flow (cwf)
* Optimal Water Flow (owf)
* Network Expansion (ne)

**Core Network Formulations**
* CNLP (convex nonlinear program used for determining network flow rates)
* MICP (relaxation-based mixed-integer convex program)
* NCNLP (non-convex nonlinear program)

## Usage at a Glance
To solve a relaxed version of the network design problem, execute the following (using the `shamir` network as an example):
```
using Ipopt
using InfrastructureModels
using JuMP
using WaterModels

ipopt = JuMP.with_optimizer(Ipopt.Optimizer, print_level=0, tol=1.0e-9, max_iter=9999)
network = WaterModels.parse_file("test/data/epanet/shamir.inp")
modifications = WaterModels.parse_file("test/data/json/shamir.json")
InfrastructureModels.update_data!(network, modifications)
result = WaterModels.run_ne(network, MICPWaterModel, ipopt, relaxed=true)
```

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
