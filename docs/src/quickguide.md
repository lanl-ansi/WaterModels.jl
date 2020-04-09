# Quick Start Guide
The following guide walks through the solution of a water network design problem using two mixed-integer linear programming formulations (MILP and MILP-R) of the network expansion (ne) problem specification.
This is to enable solution using the readily-available open-source mixed-integer linear programming solver [Cbc](https://github.com/JuliaOpt/Cbc.jl).
Other formulations rely on the availability of mixed-integer nonlinear programming solvers that support [user-defined nonlinear functions in JuMP](http://www.juliaopt.org/JuMP.jl/dev/nlp/#User-defined-Functions-1).
However, these solvers (e.g., [Juniper](https://github.com/lanl-ansi/Juniper.jl), [KNITRO](https://github.com/JuliaOpt/KNITRO.jl)) either require additional effort to register user-defined functions or are proprietary and require a commercial license.

## Installation
The latest stable release of WaterModels can be installed using the Julia package manager with
```julia
] add WaterModels
```

For the current development version, install the package using
```julia
] add WaterModels#master
```

Next, test that the package works by executing
```julia
] test WaterModels
```

Install [Cbc](https://github.com/JuliaOpt/Cbc.jl) using
```julia
] add Cbc
```

Finally, install [InfrastructureModels](https://github.com/lanl-ansi/InfrastructureModels.jl) using
```julia
] add InfrastructureModels
```

## Solving a Network Expansion (or Design) Problem
Once the above dependencies have been installed, obtain the files [`shamir.inp`](https://raw.githubusercontent.com/lanl-ansi/WaterModels.jl/master/test/data/epanet/shamir.inp) and [`shamir.json`](https://raw.githubusercontent.com/lanl-ansi/WaterModels.jl/master/test/data/json/shamir.json).
Here, `shamir.inp` is an EPANET file describing a simple seven-node, eight-link water distribution network with one reservoir, six junctions, and eight pipes.
In accord, `shamir.json` is a JSON file specifying possible pipe diameters and associated costs per unit length, per diameter setting.
The combination of data from these two files provides the required information to set up a corresponding network design problem, where the goal is to select the most cost-efficient pipe diameters while satisfying all demand in the network.

First, the diameter and cost data from the `shamir.json` problem specification must be appended to the initial data provided by the EPANET `shamir.inp` file.
To read in the EPANET data, execute the following:

```julia
using WaterModels
data = parse_file("shamir.inp")
```

Next, to read in the JSON diameter and cost data, execute the following:
```julia
modifications = parse_file("shamir.json")
```

To merge these data together, InfrastructureModels is used to update the original network data using
```julia
import InfrastructureModels
InfrastructureModels.update_data!(data, modifications)
```

Finally, the MILP formulation for the network expansion specification can be solved using
```julia
import JuMP
import Cbc

cbc = JuMP.with_optimizer(Cbc.Optimizer)
solve_ne(data, MILPWaterModel, cbc)
```

By default, only two breakpoints are used for the linear approximation of each head loss curve.
These approximations can be more finely discretized by using additional arguments to the `solve_ne` function.
For example, to employ five breakpoints per head loss curve in this formulation, the following can be executed:
```julia
solve_ne(data, MILPWaterModel, cbc, ext=Dict(:num_breakpoints=>5))
```
Note that this takes much longer to solve due to the use of more binary variables.
However, because of the finer discretization, a lower objective (design cost) can be obtained.

Instead of linear approximation, head loss curves can also be linearly outer-approximated via the MILP-R formulation.
This formulation employs less strict requirements and avoids the use of binary variables, but solutions (e.g., diameters) may not necessarily be feasible with respect to the full (nonconvex) water network physics.
To employ five outer-approximation points per (positive or negative) head loss curve in this formulation, the following can be executed
```julia
solve_ne(data, MILPRWaterModel, cbc, ext=Dict(:num_breakpoints=>5))
```

## Obtaining Results
The `run` commands in WaterModels return detailed results data in the form of a Julia `Dict`.
This dictionary can be saved for further processing as follows:
```julia
result = solve_ne(data, MILPRWaterModel, cbc, ext=Dict(:num_breakpoints=>5))
```

For example, the algorithm's runtime and final objective value can be accessed with,
```
result["solve_time"]
result["objective"]
```

The `"solution"` field contains detailed information about the solution produced by the `run` method.
For example, the following dictionary comprehension can be used to inspect the flows in the solution:
```
Dict(name => data["q"] for (name, data) in result["solution"]["pipe"])
```
Or, to obtain resistances (influenced by the selection of diameters):
```
Dict(name => data["r"] for (name, data) in result["solution"]["pipe"])
```

For more information about WaterModels result data see the [WaterModels Result Data Format](@ref) section.

## Accessing Different Formulations
The MILP formulations discussed above assume access to a mixed-integer programming (MIP) solver.
Mixed-integer nonconvex formulations can be solved with dedicated solvers, as well.
For example, the full mixed-integer nonconvex formulation for network expansion (NCNLP) can be solved via
```julia
import KNITRO

knitro = JuMP.with_optimizer(KNITRO.Optimizer)
solve_ne(data, NCNLPWaterModel, knitro)
```
and the mixed-integer convex formulation (MICP) can be solved via

```julia
solve_ne(data, MICPWaterModel, knitro)
```

## Modifying Network Data
The following example demonstrates one way to perform multiple WaterModels solves while modifing network data in Julia.
```julia
solve_ne(data, MILPRWaterModel, cbc, ext=Dict(:num_breakpoints=>5))

data["junction"]["3"]["demand"] *= 2.0
data["junction"]["4"]["demand"] *= 2.0
data["junction"]["5"]["demand"] *= 2.0

solve_ne(data, MILPRWaterModel, cbc, ext=Dict(:num_breakpoints=>5))
```
Note that the greater demands in the second problem result in an overall larger network cost.
For additional details about the network data, see the [WaterModels Network Data Format](@ref) section.

## Alternative Methods for Building and Solving Models
The following example demonstrates how to break a `solve_ne` call into separate model building and solving steps.
This allows inspection of the JuMP model created by WaterModels for the problem.
```julia
wm = build_model(data, MILPRWaterModel, WaterModels.post_ne)

print(wm.model)

result = optimize_model!(wm, cbc)
```
