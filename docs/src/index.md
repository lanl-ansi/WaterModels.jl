# WaterModels.jl Documentation

```@meta
CurrentModule = WaterModels
```

## Overview
WaterModels.jl is a Julia/JuMP package for steady state water network optimization.
It is designed to enable computational evaluation of historical and emerging water network formulations and algorithms using a common platform.
The code is engineered to decouple [Problem Specifications](@ref) (e.g., water flow, optimal water flow, network expansion) from [Network Formulations](@ref) (e.g., mixed-integer linear, mixed-integer nonlinear).
This decoupling enables the definition of a wide variety of water network optimization formulations and their comparison on common problem specifications.

## Installation
The latest stable release of WaterModels can be installed using the Julia package manager with
```julia
] add WaterModels
```

For the current development version, install the package using
```julia
] add WaterModels#master
```

Test that the package works by executing
```julia
] test WaterModels
```

Note that the WaterModels tests are comprehensive and can take as long as twenty minutes to complete.

## Usage at a Glance
At least one optimization solver is required to run WaterModels.
The solver selected typically depends on the type of problem formulation being employed.
As an example, to solve a mixed-integer linear programming (MILP) formulation of the water flow feasibility problem, the open-source mixed-integer programming solver [CBC](https://github.com/coin-or/Cbc) can be used.
Installation of the JuMP interface to CBC can be performed via the Julia package manner, i.e.,

```julia
] add Cbc
```

Then, as one example, an approximation of water flow physics for the well-known [shamir network](https://github.com/lanl-ansi/WaterModels.jl/blob/master/test/data/epanet/shamir.inp), using ten breakpoints to model each potential (or head) loss curve, can be obtained by executing

```julia
using Cbc
using WaterModels
ext = Dict(:num_breakpoints=>10)
result = run_wf("shamir.inp", MILPWaterModel, Cbc.Optimizer, ext=ext)
```
