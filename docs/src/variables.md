# Variables
The following methods provide a compositional approach for defining common variables used in water network optimization models.

## Formulation-agnostic Variables

These methods are always defined over `AbstractWaterModel`.

### Nodal Variables
```@docs
variable_head
variable_demand_flow
variable_reservoir_flow
variable_tank_flow
```

### Link Variables
```@docs
variable_des_pipe_indicator
variable_pump_head_gain
variable_pump_power
variable_pump_indicator
variable_pump_switch_off
variable_pump_switch_on
variable_regulator_indicator
variable_valve_indicator
```

## Flow-related Variables
In most of the implemented formulations, we model flow-related quantities in different ways.
The formulation-specific functions used for instantiating the associated variables are described below.
```@docs
variable_flow
```