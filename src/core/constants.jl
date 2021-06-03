"Enumerated type specifying the direction of flow along an edge."
@enum FLOW_DIRECTION FLOW_DIRECTION_POSITIVE=1 FLOW_DIRECTION_NEGATIVE=-1 FLOW_DIRECTION_UNKNOWN=0

"Ensures that JSON serialization of `FLOW_DIRECTION` returns an integer."
JSON.lower(x::FLOW_DIRECTION) = Int(x)

"Enumerated type specifying the form of the pump head curve."
@enum PUMP PUMP_QUADRATIC=0 PUMP_BEST_EFFICIENCY_POINT=1 PUMP_EPANET=2 PUMP_LINEAR_POWER=3

"Ensures that JSON serialization of `PUMP` returns an integer."
JSON.lower(x::PUMP) = Int(x)

"Enumerated type specifying the status of the component."
@enum STATUS STATUS_UNKNOWN=-1 STATUS_INACTIVE=0 STATUS_ACTIVE=1

"Ensures that JSON serialization of `STATUS` returns an integer."
JSON.lower(x::STATUS) = Int(x)

"Defines a constant for determining whether flow along a pump is considered appreciable."
const _FLOW_MIN = 6.31465679e-6

"Defines the constant for gravitational acceleration in meters per second squared"
const _GRAVITY = 9.80665

"Defines the constant for the density of water in kilograms per cubic meter"
const _DENSITY = 1000.0