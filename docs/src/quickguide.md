# Quick Start Guide
The following guide walks through the solution of a water network design (`des`) problem using two mixed-integer linear programming (MILP) formulations (PWLRD and LRD) of the problem specification.
This is to enable solution using the readily-available open-source MILP solver [HiGHS](https://github.com/jump-dev/HiGHS.jl).
Other formulations rely on the availability of mixed-integer nonlinear programming (MINLP) solvers that support [user-defined nonlinear functions in JuMP](http://www.juliaopt.org/JuMP.jl/dev/nlp/#User-defined-Functions-1).
However, these solvers (e.g., [Juniper](https://github.com/lanl-ansi/Juniper.jl), [KNITRO](https://github.com/JuliaOpt/KNITRO.jl)) either require additional effort to register user-defined functions or are proprietary and require a commercial license.

## Installation of WaterModels
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

## Installation of an Optimization Solver
At least one optimization solver is required to run WaterModels.
The solver selected typically depends on the type of problem formulation being employed. Because in this example, we will be studying the linearization-based PWLRD and LRD formulations, we will leverage the open-source MILP solver HiGHS.
Installation of the JuMP interface to HiGHS can be performed via the Julia package manager, i.e.,
```julia
] add HiGHS
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

Since we are using a linearization-based formulation of the problem, it is important to specify the partitioning of flows that parameterize the formulation.
Here, we initialize linearization flow partitions that assume a head loss error tolerance of fifty meters, with widths between flow points no greater than 1.0e-4 cubic meters per second:
```julia
set_flow_partitions_si!(data, 50.0, 1.0e-4)
```

Finally, the PWLRD formulation for the network design specification can be solved using
```julia
import HiGHS
solve_des(data, PWLRDWaterModel, HiGHS.Optimizer)
```

The above flow partitioning, however, is somewhat coarse, and the number of points in each partition is typically three, e.g.,
```julia
data["des_pipe"]["3"]["flow_partition"]
```

The relaxation can be more finely discretized by using a smaller head loss error tolerance, e.g.,
```julia
set_flow_partitions_si!(data, 5.0, 1.0e-4)
```

We can then solve the problem with the updated partitioning scheme via
```julia
import JuMP
highs = JuMP.optimizer_with_attributes(HiGHS.Optimizer, "time_limit" => 30.0)
solve_des(data, PWLRDWaterModel, highs)
```

Note that this formulation takes much longer to solve to global optimality due to the use of more binary variables.
However, because of the finer discretization, a better approximation of the physics is attained.

Instead of using piecewise-linear envelopes, head loss curves can also be simply outer-approximated via the LRD formulation.
This formulation employs less strict requirements and avoids the use of binary variables for piecewise approximation, but solutions (e.g., diameters) may not be as close to feasibility with respect to the full (nonconvex) water network physics.
To solve an LRD formulation of the problem using an even finer flow partitioning scheme (but without piecewise inner head loss approximations), the following can be executed:
```julia
set_flow_partitions_si!(data, 0.5, 1.0e-4)
solve_des(data, LRDWaterModel, HiGHS.Optimizer)
```

This relaxation of the problem turns out to converge to the known globally optimal objective value.

## Obtaining Results
For the rest of this tutorial, we will first assume a coarser relaxation by resetting the flow partitions as
```julia
set_flow_partitions_si!(data, 50.0, 1.0e-4)
```

The `solve` commands in WaterModels return detailed results data in the form of a Julia `Dict`.
This dictionary can be saved for further processing as follows:
```julia
result = solve_des(data, LRDWaterModel, HiGHS.Optimizer)
```

For example, the algorithm's runtime and final objective value can be accessed with
```julia
result["solve_time"] # Total solve time required (seconds).
result["objective"] # Final objective value (in units of the objective).
```

The `"solution"` field contains detailed information about the solution produced by the `solve` method.
For example, the following dictionary comprehension can be used to inspect the flows in the solution:
```julia
flows = Dict(name => data["q"] for (name, data) in result["solution"]["des_pipe"])
```

To determine the design pipes that were selected via the optimization, the following can be used:
```julia
pipes_selected = filter(x -> x.second["status"] == 1, result["solution"]["des_pipe"])
```

To retrieve the subset of the original pipe dataset, the following can be used:
```julia
pipes_subset = filter(x -> x.first in keys(pipes_selected), data["des_pipe"])
```

For more information about WaterModels result data see the [WaterModels Result Data Format](@ref) section.

## Accessing Different Formulations
The MILP formulations discussed above assume access to a MILP solver.
Nonconvex MINLP formulations can be solved with dedicated solvers, as well.
For example, the nonconvex MINLP formulation for design (NC) can be solved via

```julia
import KNITRO
solve_des(data, NCWaterModel, KNITRO.Optimizer)
```

## Modifying Network Data
The following example demonstrates one way to perform multiple WaterModels solves while modifying network data:
```julia
solve_des(data, LRDWaterModel, HiGHS.Optimizer)

data["demand"]["3"]["flow_min"] *= 0.5
data["demand"]["3"]["flow_max"] *= 0.5
data["demand"]["3"]["flow_nominal"] *= 0.5

solve_des(data, LRDWaterModel, HiGHS.Optimizer)
```

Note that the smaller demands in the second problem result in an overall smaller design cost.
For additional details about the network data, see the [WaterModels Network Data Format](@ref) section.

## Alternative Methods for Building and Solving Models
The following example demonstrates how to break a `solve_des` call into separate model building and solving steps.
This allows inspection of the JuMP model created by WaterModels for the problem.
```julia
wm = instantiate_model(data, LRDWaterModel, WaterModels.build_des);

println(wm.model)

result = optimize_model!(wm, optimizer = HiGHS.Optimizer)
```

## Solution Unit Conversion
The default behavior of WaterModels produces solution results in non-dimensionalized units.
To recover solutions in SI units, the following function can be used:
```julia
make_si_units!(result["solution"])
```