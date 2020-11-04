# Quick Start Guide
The following guide walks through the solution of a water network design (`des`) problem using two mixed-integer linear programming formulations (LA and LRD) of the problem specification.
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

Finally, test that the package works as expected by executing
```julia
] test WaterModels
```

## Solving a Network Design Problem
Once the above dependencies have been installed, obtain the files [`shamir.inp`](https://raw.githubusercontent.com/lanl-ansi/WaterModels.jl/master/examples/data/epanet/shamir.inp) and [`shamir.json`](https://raw.githubusercontent.com/lanl-ansi/WaterModels.jl/master/examples/data/json/shamir.json).
Here, `shamir.inp` is an EPANET file describing a simple seven-node, eight-link water distribution network with one reservoir, six demands, and eight pipes.
In accord, `shamir.json` is a JSON file specifying possible pipe diameters and associated costs per unit length, per diameter setting.
The combination of data from these two files provides the required information to set up a corresponding network design problem, where the goal is to select the most cost-efficient pipe diameters while satisfying all demand in the network.

First, the diameter and cost data from the `shamir.json` problem specification must be appended to the initial data provided by the EPANET `shamir.inp` file.
To read in the EPANET data, execute the following:

```julia
using WaterModels
data = parse_file("examples/data/epanet/shamir.inp")
```

Next, to read in the JSON diameter and cost data, execute the following:
```julia
modifications = parse_file("examples/data/json/shamir.json")
```

To merge these data together, InfrastructureModels is used to update the original network data using
```julia
WaterModels._IM.update_data!(data, modifications)
```

Finally, the LA formulation for the network design specification can be solved using
```julia
import Cbc
run_des(data, LAWaterModel, Cbc.Optimizer)
```

By default, no breakpoints are used for the linear approximation of each head loss curve, and the problem is infeasible.
These approximations can be more finely discretized by using additional arguments to the `run_des` function.
For example, to employ five breakpoints per head loss curve in this formulation, the following can be executed:
```julia
run_des(data, LAWaterModel, Cbc.Optimizer, ext=Dict(:pipe_breakpoints=>5))
```
Note that this takes much longer to solve due to the use of more binary variables.
However, because of the finer discretization, a better approximation of the physics is attained.

Instead of linear approximation, head loss curves can also be linearly outer-approximated via the LRD formulation.
This formulation employs less strict requirements and avoids the use of binary variables, but solutions (e.g., diameters) may not necessarily be feasible with respect to the full (nonconvex) water network physics.
To employ five outer approximation points per (positive or negative) head loss curve in this formulation, the following can be executed:
```julia
run_des(data, LRDWaterModel, Cbc.Optimizer, ext=Dict(:pipe_breakpoints=>5))
```

## Obtaining Results
The `run` commands in WaterModels return detailed results data in the form of a Julia `Dict`.
This dictionary can be saved for further processing as follows:
```julia
result = run_des(data, LRDWaterModel, Cbc.Optimizer, ext=Dict(:pipe_breakpoints=>5))
```

For example, the algorithm's runtime and final objective value can be accessed with,
```
result["solve_time"]
result["objective"]
```

The `"solution"` field contains detailed information about the solution produced by the `run` method.
For example, the following dictionary comprehension can be used to inspect the flows in the solution:
```
Dict(name => data["q"] for (name, data) in result["solution"]["des_pipe"])
```

For more information about WaterModels result data see the [WaterModels Result Data Format](@ref) section.

## Accessing Different Formulations
The MILP formulations discussed above assume access to a mixed-integer programming (MIP) solver.
Mixed-integer nonconvex formulations can be solved with dedicated solvers, as well.
For example, the full mixed-integer nonconvex formulation for design (NC) can be solved via
```julia
import KNITRO
run_des(data, NCWaterModel, KNITRO.Optimizer)
```

## Modifying Network Data
The following example demonstrates one way to perform multiple WaterModels solves while modifying network data:
```julia
run_des(data, LRDWaterModel, Cbc.Optimizer, ext=Dict(:pipe_breakpoints=>5))

data["demand"]["3"]["flow_rate"] *= 2.0
data["demand"]["4"]["flow_rate"] *= 2.0
data["demand"]["5"]["flow_rate"] *= 2.0

run_des(data, LRDWaterModel, Cbc.Optimizer, ext=Dict(:pipe_breakpoints=>5))
```
Note that the greater demands in the second problem result in an overall larger network cost.
For additional details about the network data, see the [WaterModels Network Data Format](@ref) section.

## Alternative Methods for Building and Solving Models
The following example demonstrates how to break a `run_des` call into separate model building and solving steps.
This allows inspection of the JuMP model created by WaterModels for the problem.
```julia
wm = instantiate_model(data, LRDWaterModel, WaterModels.build_des)

print(wm.model)

result = WaterModels._IM.optimize_model!(wm, optimizer=Cbc.Optimizer)
```
