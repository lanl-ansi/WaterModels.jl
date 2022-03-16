# WaterModels.jl Documentation

```@meta
CurrentModule = WaterModels
```

## Overview
WaterModels.jl is a Julia/JuMP package for steady-state potable water distribution network optimization.
It is designed to enable the computational evaluation of historical and emerging water network optimization formulations and algorithms using a common platform.
The code is specifically engineered to decouple [Problem Specifications](@ref) (e.g., water flow, optimal water flow, network design) from [Network Formulations](@ref) (e.g., mixed-integer linear, mixed-integer nonlinear).
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
As an example, to solve a mixed-integer linear programming (MILP) formulation of the feasible water flow (`wf`) problem, the open-source MILP solver [HiGHS](https://github.com/jump-dev/HiGHS.jl) can be used.
Installation of the JuMP interface to HiGHS can be performed via the Julia package manager, i.e.,

```julia
] add HiGHS
```

Then, as one example, a piecewise-linear, relaxation-based convexification of the physics for the well-known [Shamir (two-loop) network](https://github.com/lanl-ansi/WaterModels.jl/blob/master/examples/data/epanet/shamir.inp), using an error tolerance of one meter to model the envelope of each pipe's Hazen-Williams head loss curve, can be solved to feasibility using

```julia
using WaterModels, HiGHS

# Parse the network data from an EPANET file.
network = parse_file("examples/data/epanet/shamir.inp")

# Set linearization partitioning points that assume a head loss error tolerance of one
# meter, with widths between flow points no greater than 1.0e-4 cubic meters per second.
set_flow_partitions_si!(network, 1.0, 1.0e-4)

# Solve the corresponding relaxation of the water flow problem.
result = solve_wf(network, PWLRDWaterModel, HiGHS.Optimizer)
```

After solving the problem, results can then be analyzed, e.g.,
```julia
# The termination status of the optimization solver.
result["termination_status"]

# The flow along pipe 4 in cubic meters per second.
result["solution"]["pipe"]["4"]["q"] * result["solution"]["base_flow"]

# The total hydraulic head at node 2 in meters.
result["solution"]["node"]["2"]["h"] * result["solution"]["base_head"]

# The pressure head at node 2 in meters.
result["solution"]["node"]["2"]["p"] * result["solution"]["base_head"]
```