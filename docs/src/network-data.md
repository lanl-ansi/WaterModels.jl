# WaterModels Network Data Format

## The Network Data Dictionary
Internally, WaterModels uses a dictionary to store network data.
The dictionary uses strings as key values so it can be serialized to JSON for algorithmic data exchange.
When used, the data is assumed to be in per-unit (non-dimenisionalized) or SI units.
The network data dictionary structure is organized as follows:

```json
{
    "name": <string>,       # Name for the network model.
    "multinetwork": <bool>  # Whether or not the network data describes a multinetwork (e.g., a time series of networks).
    "per_unit": <bool>,     # Whether the data is in per-unit (non-dimensionalized) or SI units.
    "base_flow": <float>,   # Base for non-dimensionalizing volumetric flow rate. SI units are cubic meters per second.
    "base_head": <float>,   # Base for non-dimensionalizing total hydraulic head head. SI units are meters.
    "base_length": <float>, # Base for non-dimensionalizing length. SI units are meters.
    "base_mass": <float>,   # Base for non-dimensionalizing mass. SI units are kilograms.
    "base_time": <float>,   # Base for non-dimensionalizing time. SI units are seconds.
    "head_loss": <string>,  # Pipe head loss form. Can be "H-W" for Hazen-Williams or "D-W" for Darcy-Weisbach.
    "time_step": <float>,   # Time step of network data. In per-unit (non-dimensionalized) or SI units.
    "viscosity": <string>,  # Viscosity of water in the system. In per-unit (non-dimensionalized) or SI units.
    "node": {
        "1": {
            "index": <int>,           # Index of the node, which should correspond to the key.
            "name": <string>,         # Name of the node, which may or may not correspond to the index.
            "source_id": <tuple>,     # Pair (e.g., ["reservoir", "1"]) describing where the node appears in the source file.
            "status": <int>,          # Status of the node (-1 = unknown, 0 = inactive, 1 = active).
            "elevation": <float>,     # Elevation of the node above some predefined datum. SI units are meters.
            "head_min": <float>,      # Minimum total hydraulic head of the node. SI units are meters.
            "head_max": <float>,      # Maximum total hydraulic head of the node. SI units are meters.
            "head_nominal": <float>,  # Nominal total hydraulic head of the node. SI units are meters.
            "coordinates": <tuple>    # Pair of coordinates (e.g., [-1.0, 1.0]) defining the graphical location of the node.
        },
        "2": {
            ...
        },
        ...
    },
    "demand": {
        "1": {
            "index": <int>,          # Index of the demand, which should correspond to the key.
            "node": <int>,           # Index of the node to which the demand is attached.
            "name": <string>,        # Name of the demand, which may or may not correspond to the index.
            "source_id": <tuple>,    # Pair (e.g., ["demand", "1"]) describing where the demand appears in the source file.
            "status": <int>,         # Status of the demand (-1 = unknown, 0 = inactive, 1 = active).
            "dispatchable": <bool>,  # Whether or not the demand can fluctuate (true) or if it is fixed (false).
            "flow_min": <float>,     # Minimum volumetric flow rate demanded. SI units are cubic meters per second.
            "flow_max": <float>,     # Maximum volumetric flow rate demanded. SI units are cubic meters per second.
            "flow_nominal": <float>  # Nominal volumetric flow rate demanded. SI units are cubic meters per second.
        },
        "2": {
            ...
        },
        ...
    },
    "reservoir": {
        "1": {
            "index": <int>,          # Index of the reservoir, which should correspond to the key.
            "node": <int>,           # Index of the node to which the reservoir is attached.
            "name": <string>,        # Name of the reservoir, which may or may not correspond to the index.
            "source_id": <tuple>,    # Pair (e.g., ["reservoir", "1"]) describing where the reservoir appears in the source file.
            "status": <int>,         # Status of the reservoir (-1 = unknown, 0 = inactive, 1 = active).
            "dispatchable": <bool>,  # Whether or not the reservoir's head can fluctuate (true) or if it is fixed (false).
            "head_nominal": <float>  # Nominal head of the reservoir's surface. SI units are meters.
        },
        "2": {
            ...
        },
        ...
    },
    "tank": {
        "1": {
            "index": <int>,          # Index of the tank, which should correspond to the key.
            "node": <int>,           # Index of the node to which the tank is attached.
            "name": <string>,        # Name of the tank, which may or may not correspond to the index.
            "source_id": <tuple>,    # Pair (e.g., ["tank", "1"]) describing where the tank appears in the source file.
            "status": <int>,         # Status of the tank (-1 = unknown, 0 = inactive, 1 = active).
            "dispatchable": <bool>,  # Whether or not the tank's head can fluctuate (true) or if it is fixed (false).
            "diameter": <float>,     # Cross-sectional diameter of the cylindrical tank. SI units are meters.
            "min_vol": <float>,      # Minimum water volume contained by the tank. SI units are cubic meters.
            "init_level": <float>,   # Initial water level of the tank. SI units are meters.
            "min_level": <float>,    # Minimum water level of the tank. SI units are meters.
            "max_level": <float>     # Maximum water level of the tank. SI units are meters.
        },
        "2": {
            ...
        },
        ...
    },
    "pipe": {
        "1": {
            "index": <int>,               # Index of the pipe, which should correspond to the key.
            "node_fr": <int>,             # Index of the "from" node to which the pipe is connected.
            "node_to": <int>,             # Index of the "to" node to which the pipe is connected.
            "name": <string>,             # Name of the pipe, which may or may not correspond to the index.
            "source_id": <tuple>,         # Pair (e.g., ["pipe", "1"]) describing where the pipe appears in the source file.
            "status": <int>,              # Status of the pipe (-1 = unknown, 0 = inactive, 1 = active).
            "length": <float>,            # Length of the pipe. SI units are meters.
            "diameter": <float>,          # Interior (inner) diameter of the pipe. SI units are meters.
            "roughness": <float>,         # Roughness of the pipe. SI units are meters if "head_loss" is "D-W" and unitless if "H-W".
            "flow_direction": <int>,      # Direction of flow through the pipe (-1 = negative, 0 = unknown, 1 = positive).
            "flow_min": <float>,          # Minimum volumetric flow rate through the pipe. SI units are cubic meters per second.
            "flow_max": <float>,          # Maximum volumetric flow rate through the pipe. SI units are cubic meters per second.
            "flow_min_forward": <float>,  # Minimum volumetric flow rate when positively-directed. SI units are cubic meters per second.
            "flow_max_reverse": <float>,  # Maximum volumetric flow rate when negatively-directed. SI units are cubic meters per second.
            "minor_loss": <float>         # Used for modeling other minor losses of the pipe. Unitless. Currently unused.
        },
        "2": {
            ...
        },
        ...
    },
    "des_pipe": {
        "1": {
            "index": <int>,               # Index of the design pipe, which should correspond to the key.
            "node_fr": <int>,             # Index of the "from" node to which the design pipe is connected.
            "node_to": <int>,             # Index of the "to" node to which the design pipe is connected.
            "name": <string>,             # Name of the design pipe, which may or may not correspond to the index.
            "source_id": <tuple>,         # Pair (e.g., ["pipe", "1"]) describing where the design pipe appears in the source file.
            "status": <int>,              # Status of the design pipe (-1 = unknown, 0 = inactive, 1 = active).
            "length": <float>,            # Length of the design pipe. SI units are meters.
            "diameter": <float>,          # Interior (inner) diameter of the design pipe. SI units are meters.
            "roughness": <float>,         # Roughness of the design pipe. SI units are meters if "head_loss" is "D-W" and unitless if "H-W".
            "flow_direction": <int>,      # Direction of flow through the design pipe (-1 = negative, 0 = unknown, 1 = positive).
            "flow_min": <float>,          # Minimum volumetric flow rate through the design pipe. SI units are cubic meters per second.
            "flow_max": <float>,          # Maximum volumetric flow rate through the design pipe. SI units are cubic meters per second.
            "flow_min_forward": <float>,  # Minimum volumetric flow rate when positively-directed. SI units are cubic meters per second.
            "flow_max_reverse": <float>,  # Maximum volumetric flow rate when negatively-directed. SI units are cubic meters per second.
            "minor_loss": <float>,        # Used for modeling other minor losses of the design pipe. Unitless. Currently unused.
            "cost": <float>               # Cost of constructing the design pipe, if selected. Standard units are of currency.
        },
        "2": {
            ...
        },
        ...
    },
    "short_pipe": {
        "1": {
            "index": <int>,               # Index of the short pipe, which should correspond to the key.
            "node_fr": <int>,             # Index of the "from" node to which the short pipe is connected.
            "node_to": <int>,             # Index of the "to" node to which the short pipe is connected.
            "name": <string>,             # Name of the short pipe, which may or may not correspond to the index.
            "source_id": <tuple>,         # Pair (e.g., ["short_pipe", "1"]) describing where the short pipe appears in the source file.
            "status": <int>,              # Status of the short pipe (-1 = unknown, 0 = inactive, 1 = active).
            "flow_direction": <int>,      # Direction of flow through the short pipe (-1 = negative, 0 = unknown, 1 = positive).
            "flow_min": <float>,          # Minimum volumetric flow rate through the short pipe. SI units are cubic meters per second.
            "flow_max": <float>,          # Maximum volumetric flow rate through the short pipe. SI units are cubic meters per second.
            "flow_min_forward": <float>,  # Minimum volumetric flow rate when positively-directed. SI units are cubic meters per second.
            "flow_max_reverse": <float>,  # Maximum volumetric flow rate when negatively-directed. SI units are cubic meters per second.
            "minor_loss": <float>         # Used for modeling other minor losses of the short pipe. Unitless. Currently unused.
        },
        "2": {
            ...
        },
        ...
    },
    "pump": {
        "1": {
            "index": <int>,                   # Index of the pump, which should correspond to the key.
            "node_fr": <int>,                 # Index of the "from" node to which the pump is connected.
            "node_to": <int>,                 # Index of the "to" node to which the pump is connected.
            "name": <string>,                 # Name of the pump, which may or may not correspond to the index.
            "source_id": <tuple>,             # Pair (e.g., ["pump", "1"]) describing where the pump appears in the source file.
            "status": <int>,                  # Status of the pump (-1 = unknown, 0 = inactive, 1 = active).
            "flow_direction": <int>,          # Direction of flow through the pump (-1 = negative, 0 = unknown, 1 = positive).
            "flow_min": <float>,              # Minimum volumetric flow rate through the pump. SI units are cubic meters per second.
            "flow_max": <float>,              # Maximum volumetric flow rate through the pump. SI units are cubic meters per second.
            "flow_min_forward": <float>,      # Minimum volumetric flow rate when positively-directed. SI units are cubic meters per second.
            "flow_max_reverse": <float>,      # Maximum volumetric flow rate when negatively-directed. SI units are cubic meters per second.
            "head_curve_form": <int>,         # Form of the head curve function used to model the pump. PUMP_QUADRATIC (0) models the head
                                              # gain as a quadratic function. PUMP_BEST_EFFICIENCY_POINT (1) models the head gain as a
                                              # quadratic function parameterized using the pump's best efficiency point. PUMP_EPANET (2)
                                              # models the head gain similarly to the model of EPANET, i.e., a + b * q^c, where a, b, and
                                              # c are fixed constants computed from flow-head gain data present in the "head_curve" list.
                                              # PUMP_LINEAR_POWER (3) models the head gain function as PUMP_QUADRATIC but power linearly.
            "head_curve": <List[...]>,        # List of tuples describing the head gain of the pump, where the first element of the tuple
                                              # corresponds to a volumetric flow rate (in per-unit or SI units) and the second element of
                                              # corresponds to the head gain resulting from that flow rate (in per-unit or SI units).
            "efficiency_curve": <List[...]>,  # List of tuples describing the efficiency of the pump, where the first element of the tuple
                                              # corresponds to a volumetric flow rate (in per-unit or SI units) and the second element
                                              # corresponds to the pump's power efficiency at that flow rate (a unitless quantity).
            "energy_price": <float>,          # Cost of consuming energy (in per-unit or, for SI units, currency per Joule).
            "power_fixed": <float>,           # When modeling the pump using PUMP_LINEAR_POWER for "head_curve_form", this quantity
                                              # represents the power consumed by the pump at zero volumetric flow. SI units are Watts.
            "power_per_unit_flow": <float>,   # When modeling the pump using PUMP_LINEAR_POWER for "head_curve_form", this quantity
                                              # represents the power consumed per additional unit of flow. SI units are Watts per
                                              # (cubic meters per second).
        },
        "2": {
            ...
        },
        ...
    },
    "valve": {
        "1": {
            "index": <int>,               # Index of the valve, which should correspond to the key.
            "node_fr": <int>,             # Index of the "from" node to which the valve is connected.
            "node_to": <int>,             # Index of the "to" node to which the valve is connected.
            "name": <string>,             # Name of the valve, which may or may not correspond to the index.
            "source_id": <tuple>,         # Pair (e.g., ["pipe", "1"]) describing where the valve appears in the source file.
            "status": <int>,              # Status of the valve (-1 = unknown, 0 = inactive, 1 = active).
            "flow_direction": <int>,      # Direction of flow through the valve (-1 = negative, 0 = unknown, 1 = positive).
            "flow_min": <float>,          # Minimum volumetric flow rate through the valve. SI units are cubic meters per second.
            "flow_max": <float>,          # Maximum volumetric flow rate through the valve. SI units are cubic meters per second.
            "flow_min_forward": <float>,  # Minimum volumetric flow rate when positively-directed. SI units are cubic meters per second.
            "flow_max_reverse": <float>,  # Maximum volumetric flow rate when negatively-directed. SI units are cubic meters per second.
            "minor_loss": <float>         # Used for modeling other minor losses of the valve. Unitless. Currently unused.
        },
        "2": {
            ...
        },
        ...
    },
    "regulator": {
        "1": {
            "index": <int>,               # Index of the regulator, which should correspond to the key.
            "node_fr": <int>,             # Index of the "from" node to which the regulator is connected.
            "node_to": <int>,             # Index of the "to" node to which the regulator is connected.
            "name": <string>,             # Name of the regulator, which may or may not correspond to the index.
            "source_id": <tuple>,         # Pair (e.g., ["regulator", "1"]) describing where the regulator appears in the source file.
            "status": <int>,              # Status of the regulator (-1 = unknown, 0 = inactive, 1 = active).
            "flow_direction": <int>,      # Direction of flow through the regulator (-1 = negative, 0 = unknown, 1 = positive).
            "diameter": <float>,          # (Possibly artificial) diameter of the regulator. SI units are meters. Currently unused.
            "setting": <float>,           # Setting of the total hydraulic head at the node downstream of the regulator. SI units are meters.
            "flow_min": <float>,          # Minimum volumetric flow rate through the regulator. SI units are cubic meters per second.
            "flow_max": <float>,          # Maximum volumetric flow rate through the regulator. SI units are cubic meters per second.
            "flow_min_forward": <float>,  # Minimum volumetric flow rate when positively-directed. SI units are cubic meters per second.
            "flow_max_reverse": <float>,  # Maximum volumetric flow rate when negatively-directed. SI units are cubic meters per second.
            "minor_loss": <float>         # Used for modeling other minor losses of the regulator. Unitless. Currently unused.
        },
        "2": {
            ...
        },
        ...
    }
}
```

## Multinetwork (Time Series) Data
[The Network Data Dictionary](@ref) outlines the hierarchy for a snapshot of data for a water network at some moment in time.
In practice, water network optimization problems often model evolution of the system over time.
To model temporal aspects, a similar data hierarchy is used.
These "multinetwork" dictionaries are of the form
```json
{
    "name": <string>,
    "multinetwork": true,
    "per_unit": <bool>,
    ...,
    "nw": {
        "1": {
            "name": <string>,
            "time_step": <float>,
            "node": {
                "1": {
                    ...
                },
                ...
            },
            ...
        },
        "2": {
            ...
        },
        ...
    }
}
```
That is, the _component_ dictionaries, i.e., `["tank", "regulator", "pump", "name", "des_pipe", "demand", "reservoir", "short_pipe", "node", "valve", "pipe"]`, as well as the `"time_step"` and `"name`" parameters are stored within an indexed `"nw"` subdictionary.
Here, the `"1"` subdictionary of `"nw"` might represent the network at the first time index in the model of interest, and `"2"` might represent the network at the second time index.