# WaterModels.jl
<img src="https://lanl-ansi.github.io/WaterModels.jl/dev/assets/logo.svg" align="left" width="200" alt="WaterModels Logo">

<a href="https://lanl-ansi.github.io/WaterModels.jl/latest/"><img align="top" src="https://img.shields.io/badge/docs-latest-blue.svg" alt="Latest Documentation Status"></a> <a href="https://github.com/lanl-ansi/WaterModels.jl/actions?query=workflow%3ACI"><img src="https://github.com/lanl-ansi/WaterModels.jl/workflows/CI/badge.svg?branch=master" align="top" alt="Development Build Status"></a> <a href="https://codecov.io/gh/lanl-ansi/WaterModels.jl"><img align="top" src="https://codecov.io/gh/lanl-ansi/WaterModels.jl/branch/master/graph/badge.svg" alt="Code Coverage Status"></a>

WaterModels.jl is a Julia/JuMP package for steady state water network optimization.
It is designed to enable computational evaluation of historical and emerging water network formulations and algorithms using a common platform.
The code is engineered to decouple problem specifications (e.g., network design, optimal water flow) from network optimization formulations (e.g., mixed-integer linear, mixed-integer nonlinear).
This decoupling enables the definition of a wide variety of optimization formulations and their comparison on common problem specifications.

**Core Problem Specifications**
* Water Flow (`wf` and `mn_wf`) - obtain feasible flows and heads for a network
* Optimal Water Flow (`owf` and `mn_owf`) - minimize the cost of network operation
* Network Design (`des`) - minimize the cost of network design
* Network Expansion (`ne`) - minimize the cost of network expansion

**Core Network Formulations**
* NC - nonconvex nonlinear formulation
* NCD - nonconvex nonlinear direction-based formulation
* CRD - continuous (convex) relaxation- and direction-based formulation
* LRD - linear relaxation- and direction-based formulation
* PWLRD - piecewise linear relaxation- and direction-based formulation
* LA - linear approximation-based formulation

## Documentation
The package [documentation](https://lanl-ansi.github.io/WaterModels.jl/latest) includes a [quick start guide](https://lanl-ansi.github.io/WaterModels.jl/latest/quickguide).
Be advised that aside from the quick start guide, documentation is under development and may currently be inaccurate.

## Development
Community-driven development and enhancement of WaterModels is welcomed and encouraged.
Please feel free to fork this repository and share your contributions to the master branch with a pull request.
That said, it is important to keep in mind the code quality requirements and scope of WaterModels before preparing a contribution.
See [CONTRIBUTING.md](https://github.com/lanl-ansi/WaterModels.jl/blob/master/CONTRIBUTING.md) for code contribution guidelines.

## Acknowledgments
This work is currently supported by the Advanced Grid Modeling Program within the U.S. Department of Energy under the project "Coordinated Planning and Operation of Water and Power Infrastructures for Increased Resilience and Reliability."
Work at Los Alamos National Laboratory is conducted under the auspices of the National Nuclear Security Administration of the U.S. Department of Energy under Contract No. 89233218CNA000001.
Previous work was supported by the Los Alamos National Laboratory Directed Research and Development program under the project ["Adaptation Science for Complex Natural-engineered Systems"](http://www.lanl.gov/projects/nesma) (20180033DR).
It is also supported by the [Advanced Network Science Initiative](https://lanl-ansi.github.io) at Los Alamos National Laboratory.

The primary developer is [Byron Tasseff](https://github.com/tasseff) with support from the following contributors:
- [Clayton Barrows](https://github.com/claytonpbarrows), National Renewable Energy Laboratory
- [Russell Bent](https://github.com/rb004f), Los Alamos National Laboratory
- [Carleton Coffrin](https://github.com/ccoffrin), Los Alamos National Laboratory
- Donatella Pasqualini, Los Alamos National Laboratory
- [Devon Sigler](https://github.com/dsigler1234), National Renewable Energy Laboratory
- [Jonathan Stickel](https://github.com/jjstickel), National Renewable Energy Laboratory

## License
This code is provided under a [modified BSD license](https://github.com/lanl-ansi/WaterModels.jl/blob/master/LICENSE.md) as part of the Multi-Infrastructure Control and Optimization Toolkit (MICOT), LA-CC-13-108.
