"Enumerated type specifying the direction of flow along an edge."
@enum FLOW_DIRECTION POSITIVE=1 NEGATIVE=-1 UNKNOWN=0

"Defines a constant for determining whether flow along a pump or valve is
considered appreciable. If not, the heads at adjacent nodes will be decoupled."
const _q_eps = 6.31465679e-6

"Models with two positive flow variables, one for each direction."
abstract type AbstractDirectedFlowModel <: AbstractWaterModel end

"Models with one flow variable, with the sign of the variable indicating the flow direction."
abstract type AbstractUndirectedFlowModel <: AbstractWaterModel end

"Models derived from AbstractDirectedFlowModel"
abstract type AbstractMICPModel <: AbstractDirectedFlowModel end
mutable struct MICPWaterModel <: AbstractMICPModel @wm_fields end
mutable struct MICPRWaterModel <: AbstractMICPModel @wm_fields end
abstract type AbstractMILPRModel <: AbstractDirectedFlowModel end
mutable struct MILPRWaterModel <: AbstractMILPRModel @wm_fields end

"Models derived from AbstractUndirectedFlowModel"
abstract type AbstractNLPModel <: AbstractUndirectedFlowModel end
mutable struct NLPWaterModel <: AbstractNLPModel @wm_fields end
abstract type AbstractMILPModel <: AbstractUndirectedFlowModel end
mutable struct MILPWaterModel <: AbstractMILPModel @wm_fields end

"Directed models that also use direction variables."
AbstractUndirectedModel = Union{AbstractNLPModel,AbstractMILPModel}
AbstractDirectedModel = Union{AbstractMICPModel,AbstractMILPRModel}
AbstractNonlinearModel = Union{AbstractNLPModel,AbstractMICPModel}
