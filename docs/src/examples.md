# WaterModels Examples

The `examples` directory contains two water network optimization instances that have been developed or modified from two literature instances.

The first is the famous "two-loop" water network design instance.
(It is sometimes titled after one of its authors as `shamir`.)
This design instance dates back to 1977, first appearing in [1].
The globally optimal design cost is known to be \$419,000.
Solutions of this instance using various formulation types and assumptions appeared in the [Quick Start Guide](@ref).
As an example, it can be solved using a linear relaxation-based formulation (`LRDWaterModel`) via the following:
```julia
using WaterModels
import HiGHS

data = parse_file("examples/data/json/shamir.json")
set_flow_partitions_si!(data, 0.5, 1.0e-4)
result = solve_des(data, LRDWaterModel, HiGHS.Optimizer)
```

The second is a modified version of the popular `van_zyl` optimal water flow instance, which first appeared in [2] and is also named after one of that article's authors.
Unlike the design problem, this problem has temporal aspects.
It can be constructed and solved (e.g., using the `LRDWaterModel` formulation) using the following:
```julia
using WaterModels
import HiGHS
import JuMP

data = parse_file("examples/data/epanet/van_zyl.inp")
data_mn = WaterModels.make_multinetwork(data)
set_flow_partitions_si!(data_mn, 1.0, 1.0e-4)
highs = JuMP.optimizer_with_attributes(HiGHS.Optimizer, "time_limit" => 60.0)
result = solve_mn_owf(data_mn, LRDWaterModel, highs)
```
Note that results are presented in an automatically-applied per-unit system.
The instance is also challenging, and only a feasible solution is returned within the time limit for the script above.

## References
[1] Alperovits, E., & Shamir, U. (1977). Design of optimal water distribution systems. _Water Resources Research_, _13_(6), 885-900.

[2] Van Zyl, J. E., Savic, D. A., & Walters, G. A. (2004). Operational optimization of water distribution systems using a hybrid genetic algorithm. _Journal of Water Resources Planning and Management_, _130_(2), 160-170.