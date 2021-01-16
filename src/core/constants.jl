"Enumerated type specifying the direction of flow along an edge."
@enum FLOW_DIRECTION POSITIVE=1 NEGATIVE=-1 UNKNOWN=0

"Ensures that JSON serialization of `FLOW_DIRECTION` returns an integer."
JSON.lower(x::FLOW_DIRECTION) = Int(x)

"Enumerated type specifying the form of the pump head curve."
@enum HEAD_CURVE_FORM QUADRATIC=0 BEST_EFFICIENCY_POINT=1 EPANET=2

"Ensures that JSON serialization of `HEAD_CURVE_FORM` returns an integer."
JSON.lower(x::HEAD_CURVE_FORM) = Int(x)

"Defines a constant for determining whether flow along a pump is considered appreciable."
const _FLOW_MIN = 6.31465679e-6

"Defines the constant for gravitational acceleration in meters per second squared"
const _GRAVITY = 9.80665

"Defines the constant for the density of water in kilograms per cubic meter"
const _DENSITY = 1000.0