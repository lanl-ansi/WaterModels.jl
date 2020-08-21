# WaterModels.jl Documentation

```@meta
CurrentModule = WaterModels
```

## Overview
WaterModels.jl is a Julia/JuMP package for steady state potable water distribution network optimization.
It is designed to enable the computational evaluation of historical and emerging water network optimization formulations and algorithms using a common platform.
The code is engineered to decouple [Problem Specifications](@ref) (e.g., water flow, optimal water flow, network design) from [Network Formulations](@ref) (e.g., mixed-integer linear, mixed-integer nonlinear).
This decoupling enables the definition of a wide variety of water network optimization formulations and their comparison across several common problem specifications.

## Installation
The latest stable release of WaterModels can be installed using the Julia package manager with
```julia
] add WaterModels
```

For the current development version, install the package using
```julia
] add WaterModels#master
```

Finally, test that the package works as expected by executing
```julia
] test WaterModels
```

## Usage at a Glance
At least one optimization solver is required to run WaterModels.
The solver selected typically depends on the type of problem formulation being employed.
As an example, to solve a mixed-integer linear programming (MILP) formulation of the feasible water flow (`wf`) problem, the open-source mixed-integer programming solver [CBC](https://github.com/coin-or/Cbc) can be used.
Installation of the JuMP interface to CBC can be performed via the Julia package manager, i.e.,

```julia
] add Cbc
```

Then, as one example, a piecewise-linear approximation of the physics for the well-known [Shamir (two-loop) network](https://github.com/lanl-ansi/WaterModels.jl/blob/master/examples/data/epanet/shamir.inp), using ten breakpoints to model each pipe's Hazen-Williams head loss curve, can be solved to feasibility using

```julia
using WaterModels, Cbc
ext = Dict(:pipe_breakpoints=>10) # Defines additional parameters used in model construction.
result = run_wf("examples/data/epanet/shamir.inp", MILPWaterModel, Cbc.Optimizer; ext=ext)
```

After solving the problem, its results can then be analyzed, e.g.,
```julia
# The termination status of the optimization solver.
result["termination_status"]

# The flow along pipe 4, in cubic meters per second.
result["solution"]["pipe"]["4"]["q"]

# The total hydraulic head at junction (node) 2, in meters.
result["solution"]["node"]["2"]["h"]

# The pressure head at junction (node) 2, in meters.
result["solution"]["node"]["2"]["p"]
```
