"Enumerated type specifying the direction of flow along an edge."
@enum FLOW_DIRECTION begin
    FLOW_DIRECTION_POSITIVE = 1 # The flow is directed from i to j.
    FLOW_DIRECTION_NEGATIVE = -1 # The flow is directed from j to i.
    FLOW_DIRECTION_UNKNOWN = 0 # The flow direction is unknown.
end

"Ensures that JSON serialization of `FLOW_DIRECTION` returns an integer."
JSON.lower(x::FLOW_DIRECTION) = Int(x)

"Enumerated type specifying the method used for modeling a pump."
@enum PUMP begin
    PUMP_QUADRATIC = 0 # Head gain takes a quadratic form.
    PUMP_BEST_EFFICIENCY_POINT = 1 # Head gain takes a quadratic best efficiency form.
    PUMP_EPANET = 2 # Head gain takes the form used by EPANET.
    PUMP_LINEAR_POWER = 3 # Power is modeled linearly and head gain quadratic.
    PUMP_CONSTANT_POWER = 4 # Power is constant when the pump is active.
end

"Ensures that JSON serialization of `PUMP` returns an integer."
JSON.lower(x::PUMP) = Int(x)

"Enumerated type specifying the status of the component."
@enum STATUS begin
    STATUS_UNKNOWN = -1 # The status of the component is unknown (i.e., on or off).
    STATUS_INACTIVE = 0 # The status of the component is inactive (i.e., off or removed).
    STATUS_ACTIVE = 1 # The status of the component is active (i.e., on or present).
end

"Ensures that JSON serialization of `STATUS` returns an integer."
JSON.lower(x::STATUS) = Int(x)

"Defines a constant for determining whether flow along a pump is considered appreciable."
const _FLOW_MIN = 6.31465679e-6 # In cubic meters per second.

"Defines the constant for gravitational acceleration."
const _GRAVITY = 9.80665 # In meters per second squared.

"Defines the constant for the density of water."
const _DENSITY = 1000.0 # In kilograms per cubic meter.

"Defines a convenient vector for node-connected component types."
const _NODE_CONNECTED_COMPONENTS = ["demand", "reservoir", "tank"]

"Defines a convenient vector for node-connecting component types."
const _LINK_COMPONENTS = ["pipe", "pump", "ne_pump", "regulator", "short_pipe", "valve", "des_pipe", "ne_short_pipe"]
