# Quick Start Guide
The following guide walks through the solution of a water network design (`des`) problem using two mixed-integer linear programming formulations (PWLRD and LRD) of the problem specification.
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
Once the above dependencies have been installed, obtain the file [`shamir.json`](https://raw.githubusercontent.com/lanl-ansi/WaterModels.jl/master/examples/data/json/shamir.json).
Here, `shamir.json` is a JSON file specifying the network, as well as possible pipe diameters and associated costs, per diameter setting.
The file provides the required information to set up a corresponding network design problem, where the goal is to select the most cost-efficient pipe diameters while satisfying all demand in the network.

To read in the data, execute the following:

```julia
using WaterModels
data = parse_file("examples/data/json/shamir.json")
```

Finally, the PWLRD formulation for the network design specification can be solved using
```julia
import Cbc
solve_des(data, PWLRDWaterModel, Cbc.Optimizer)
```

By default, two breakpoints are used for the linear approximation of each directed head loss curve.
These approximations can be more finely discretized by using additional arguments to the `solve_des` function.
For example, to employ five breakpoints per head loss curve in this formulation, the following can be executed:
```julia
import JuMP
cbc = JuMP.optimizer_with_attributes(Cbc.Optimizer, "seconds" => 30.0)
solve_des(data, PWLRDWaterModel, cbc, ext=Dict(:pipe_breakpoints=>5))
```
Note that this formulation takes much longer to solve to global optimality due to the use of more binary variables.
However, because of the finer discretization, a better approximation of the physics is attained.

Instead of using piecewise-linear envelopes, head loss curves can also be simply outer-approximated via the LRD formulation.
This formulation employs less strict requirements and avoids the use of binary variables for piecewise approximation, but solutions (e.g., diameters) may not be as close to feasibility with respect to the full (nonconvex) water network physics.
To employ five outer approximation points per directed head loss curve in this formulation, the following can be executed:
```julia
solve_des(data, LRDWaterModel, Cbc.Optimizer, ext=Dict(:pipe_breakpoints=>5))
```
This relaxation of the problem turns out to converge to the known globally optimal objective value.

## Obtaining Results
The `run` commands in WaterModels return detailed results data in the form of a Julia `Dict`.
This dictionary can be saved for further processing as follows:
```julia
result = solve_des(data, LRDWaterModel, Cbc.Optimizer)
```

For example, the algorithm's runtime and final objective value can be accessed with,
```
result["solve_time"]
result["objective"]
```

The `"solution"` field contains detailed information about the solution produced by the `run` method.
For example, the following dictionary comprehension can be used to inspect the flows in the solution:
```
flows = Dict(name => data["q"] for (name, data) in result["solution"]["des_pipe"])
```

To determine the design pipes that were selected via the optimization, the following can be used:
```
pipes_selected = filter(x -> x.second["status"] == 1, result["solution"]["des_pipe"])
```

To retrieve the subset of the original pipe dataset, the following can be used:
```
pipes_subset = filter(x -> x.first in keys(pipes_selected), data["des_pipe"])
```

For more information about WaterModels result data see the [WaterModels Result Data Format](@ref) section.

## Accessing Different Formulations
The MILP formulations discussed above assume access to a mixed-integer programming (MIP) solver.
Mixed-integer nonconvex formulations can be solved with dedicated solvers, as well.
For example, the full mixed-integer nonconvex formulation for design (NC) can be solved via
```julia
import KNITRO
solve_des(data, NCWaterModel, KNITRO.Optimizer)
```

## Modifying Network Data
The following example demonstrates one way to perform multiple WaterModels solves while modifying network data:
```julia
solve_des(data, LRDWaterModel, Cbc.Optimizer, ext=Dict(:pipe_breakpoints=>3))

data["demand"]["3"]["flow_min"] *= 0.5
data["demand"]["3"]["flow_max"] *= 0.5
data["demand"]["3"]["flow_nominal"] *= 0.5

solve_des(data, LRDWaterModel, Cbc.Optimizer, ext=Dict(:pipe_breakpoints=>3))
```
Note that the smaller demands in the second problem result in an overall smaller network cost.
For additional details about the network data, see the [WaterModels Network Data Format](@ref) section.

## Alternative Methods for Building and Solving Models
The following example demonstrates how to break a `solve_des` call into separate model building and solving steps.
This allows inspection of the JuMP model created by WaterModels for the problem.
```julia
wm = instantiate_model(data, LRDWaterModel, WaterModels.build_des);

print(wm.model)

result = WaterModels._IM.optimize_model!(wm, optimizer=Cbc.Optimizer)
```