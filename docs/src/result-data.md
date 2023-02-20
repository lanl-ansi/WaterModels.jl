# WaterModels Result Data Format

## The Result Data Dictionary
WaterModels uses a dictionary to organize the results of a `solve_` command.
The dictionary uses strings as key values so it can be serialized to JSON for algorithmic data exchange.
The data dictionary organization is designed to be consistent with [The Network Data Dictionary](@ref).

At the top level, the results data dictionary is structured as follows:
```json
{
    "optimizer": <string>,        # name of the solver used to solve the model
    "termination_status": <type>, # optimizer status at termination
    "dual_status": <type>,        # optimizer dual status at termination
    "primal_status": <type>,      # optimizer primal status at termination
    "solve_time": <float>,        # reported solve time (in seconds)
    "objective": <float>,         # the final evaluation of the objective function
    "objective_lb": <float>,      # the final lower bound of the objective function (if available)
    "solution": {...}             # complete solution information (details below)
}
```

## Solution Data
The `"solution"` subdictionary provides detailed information about the problem solution produced by the `solve_` command.
The solution is organized similarly to [The Network Data Dictionary](@ref) with the same nested structure and parameter names, when available.
For example, for a single-network problem, `result["solution"]["pipe"]["1"]` reports all the solution values associated with the pipe at index `"1"`, e.g.,
```json
{
    "qn": 0.0,
    "qp": 2.0,
    "dhn": 0.0,
    "q": 2.0,
    "dhp": 0.033770,
    "y": 1.0
}
```
Similarly, for a multinetwork problem, as an example, `result["nw"]["3"]["pump"]["2"]` reports the solution values associated with the pump at pump index `"2"` and time index `"3"`, e.g.,
```json
{
    "qn": 0.0,
    "c": 0.16850289352749692,
    "g": 0.8372647243015188,
    "P": 0.00016374541106904265,
    "status": 1.0,
    "qp": 0.36694422549881334,
    "q": 0.36694422549881334,
    "E": 1.8023959936938814e-6,
    "y": 1.0
}
```

### Solution Data Schema
By default, all solution data are reported in per-unit (non-dimensionalized) units.
Below are common outputs of the implemented optimization models, which are sometimes based on flow directionality:
```json
{
    "multiinfrastructure": <bool>  # Whether or not the solution data is part of a broader multi-infrastructure solution.
    "multinetwork": <bool>         # Whether or not the network data describes a multinetwork (e.g., a time series of networks).
    "per_unit": <bool>,            # Whether the data is in per-unit (non-dimensionalized) or SI units.
    "base_flow": <float>,          # Base for non-dimensionalizing volumetric flow rate. SI units are cubic meters per second.
    "base_head": <float>,          # Base for non-dimensionalizing total hydraulic head head. SI units are meters.
    "base_length": <float>,        # Base for non-dimensionalizing length. SI units are meters.
    "base_mass": <float>,          # Base for non-dimensionalizing mass. SI units are kilograms.
    "base_time": <float>,          # Base for non-dimensionalizing time. SI units are seconds.
    "node": {
        "1": {
            "h": <float>,  # Total hydraulic head of the node. SI units are meters.
            "p": <float>   # Pressure head of the node. SI units are meters.
        },
        "2": {
            ...
        },
        ...
    },
    "demand": {
        "1": {
            "q": <float>  # Demanded volumetric flow rate at the demand point. SI units are cubic meters per second.
        },
        "2": {
            ...
        },
        ...
    },
    "reservoir": {
        "1": {
            "q": <float>  # Outgoing volumetric flow rate from the reservoir. SI units are cubic meters per second.
        },
        "2": {
            ...
        },
        ...
    },
    "tank": {
        "1": {
            "V": <float>, # Volume of water contained by the tank. SI units are cubic meters.
            "q": <float>  # Outgoing volumetric flow rate from the reservoir. SI units are cubic meters per second.
        },
        "2": {
            ...
        },
        ...
    },
    "pipe": {
        "1": {
            "q": <float>,    # Volumetric flow rate transported through the pipe. SI units are cubic meters per second.
            "qp": <float>,   # Volumetric flow rate transported in the positive direction. SI units are cubic meters per second.
            "qn": <float>,   # Volumetric flow rate transported in the negative direction. SI units are cubic meters per second.
            "dhp": <float>,  # Total hydraulic head decrease in the positive direction of flow. SI units are meters.
            "dhn": <float>,  # Total hydraulic head decrease in the negative direction of flow. SI units are meters.
            "y": <float>     # Flow direction, i.e., one if flow is transported _from_ "node_fr" and zero otherwise.
        },
        "2": {
            ...
        },
        ...
    },
    "des_pipe": {
        "1": {
            "q": <float>,      # Volumetric flow rate transported through the design pipe. SI units are cubic meters per second.
            "qp": <float>,     # Volumetric flow rate transported in the positive direction. SI units are cubic meters per second.
            "qn": <float>,     # Volumetric flow rate transported in the negative direction. SI units are cubic meters per second.
            "dhp": <float>,    # Total hydraulic head decrease in the positive direction of flow. SI units are meters.
            "dhn": <float>,    # Total hydraulic head decrease in the negative direction of flow. SI units are meters.
            "y": <float>,      # Flow direction, i.e., one if flow is transported _from_ "node_fr" and zero otherwise.
            "status": <float>  # Pipe construction status, i.e., one if the pipe is constructed and zero otherwise.
        },
        "2": {
            ...
        },
        ...
    },
    "short_pipe": {
        "1": {
            "q": <float>,   # Volumetric flow rate transported through the short pipe. SI units are cubic meters per second.
            "qp": <float>,  # Volumetric flow rate transported in the positive direction. SI units are cubic meters per second.
            "qn": <float>,  # Volumetric flow rate transported in the negative direction. SI units are cubic meters per second.
            "y": <float>    # Flow direction, i.e., one if flow is transported _from_ "node_fr" and zero otherwise.
        },
        "2": {
            ...
        },
        ...
    },
    "pump": {
        "1": {
            "q": <float>,      # Volumetric flow rate transported through the pump. SI units are cubic meters per second.
            "qp": <float>,     # Volumetric flow rate transported in the positive direction. SI units are cubic meters per second.
            "qn": <float>,     # Volumetric flow rate transported in the negative direction. SI units are cubic meters per second.
            "y": <float>,      # Flow direction, i.e., one if flow is transported _from_ "node_fr" and zero otherwise.
            "g": <float>,      # Head gain (increase) from "node_fr" to "node_to" resulting from the pump. SI units are meters.
            "c": <float>,      # Cost of operating the pump over the time index (step) of interest. Standard units are currency.
            "P": <float>,      # Power consumed by the pump over the time index (step) of interest. SI units are Watts.
            "E": <float>,      # Energy consumption of the pump over the time index (step) of interest. SI units are Joules.
            "status": <float>  # Status of the pump, i.e., one if the pump is active and zero otherwise.
        },
        "2": {
            ...
        },
        ...
    },
    "valve": {
        "1": {
            "q": <float>,      # Volumetric flow rate transported through the valve. SI units are cubic meters per second.
            "qp": <float>,     # Volumetric flow rate transported in the positive direction. SI units are cubic meters per second.
            "qn": <float>,     # Volumetric flow rate transported in the negative direction. SI units are cubic meters per second.
            "y": <float>       # Flow direction, i.e., one if flow is transported _from_ "node_fr" and zero otherwise.
            "status": <float>  # Status of the valve, i.e., one if the valve is opened and zero otherwise.
        },
        "2": {
            ...
        },
        ...
    },
    "regulator": {
        "1": {
            "q": <float>,      # Volumetric flow rate transported through the regulator. SI units are cubic meters per second.
            "qp": <float>,     # Volumetric flow rate transported in the positive direction. SI units are cubic meters per second.
            "qn": <float>,     # Volumetric flow rate transported in the negative direction. SI units are cubic meters per second.
            "y": <float>       # Flow direction, i.e., one if flow is transported _from_ "node_fr" and zero otherwise.
            "status": <float>  # Status of the regulator, i.e., one if the regulator is active and zero otherwise.
        },
        "2": {
            ...
        },
        ...
    }
}
```


## Transforming Solution Data
Because the data dictionary and the solution dictionary have the same structure, the InfrastructureModels `update_data!` helper function can be used to update a data dictionary with values from a solution, e.g.,
```julia
import InfrastructureModels

InfrastructureModels.update_data!(
  data["nw"]["3"]["pump"]["1"],
  result["solution"]["nw"]["3"]["pump"]["1"]
)
```

Note that, by default, all results are reported in a per-unit (non-dimensionalized) system.
Additional data from WaterModels can be used to convert such data back to their dimensionalized forms.
For example, the code block below translates a per-unit pump flow rate to SI units, then the more conventional units of liters per second.
```julia
# Get a pump volumetric flow rate solution in the per-unit system.
flow_per_unit = result["solution"]["nw"]["3"]["pump"]["1"]["q"]

# Get the per-unit scalar used to convert back to SI units.
base_flow = data["base_flow"]

# Compute the volumetric flow rate in SI units (cubic meters per second).
base_flow * flow_per_unit

# Compute the volumetric flow rate in liters per second.
base_flow * flow_per_unit * 1000.0
```

Note also that per-unit base quantities are also reported in the solution data dictionary for convenience, e.g.,
```julia
# Should return `true`.
result["solution"]["base_flow"] == data["base_flow"]
```
For convenience, solution data can be transformed to SI units using
```julia
make_si_units!(result["solution"])
```