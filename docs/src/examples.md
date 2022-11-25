# WaterModels Examples

The `examples` directory contains two water network optimization instances that have been developed or modified from two literature instances.
The first is the famous "two-loop" water network design instance.
(It is sometimes titled after one of its authors as `shamir`.)
This design instance dates back to 1977, first appearing in [1].
The globally optimal design cost is known to be $419,000.
Solutions of this instance using various formulation types and assumptions appeared in [Quick Start Guide](@ref).
As an example, it can be solved using a linear relaxation-based formulation (`LRD`) via the following:
```julia
using WaterModels
import HiGHS

data = parse_file("examples/data/json/shamir.json")
set_flow_partitions_si!(data, 0.5, 1.0e-4)
result = solve_des(data, LRDWaterModel, HiGHS.Optimizer)
```

## References
[1] Alperovits, E., & Shamir, U. (1977). Design of optimal water distribution systems. _Water Resources Research_, _13_(6), 885-900.