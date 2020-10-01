"Enumerated type specifying the direction of flow along an edge."
@enum FLOW_DIRECTION POSITIVE=1 NEGATIVE=-1 UNKNOWN=0

"Defines a constant for determining whether flow along a pump or valve is
considered appreciable. If not, the heads at adjacent nodes will be decoupled."
const _q_eps = 6.31465679e-6

"Models derived from AbstractDirectedFlowModel"
mutable struct CRDWaterModel <: AbstractWaterModel @wm_fields end
mutable struct LRDWaterModel <: AbstractWaterModel @wm_fields end
mutable struct QRDWaterModel <: AbstractWaterModel @wm_fields end
mutable struct CQRDWaterModel <: AbstractWaterModel @wm_fields end

"Models derived from AbstractUndirectedFlowModel"
mutable struct NCWaterModel <: AbstractWaterModel @wm_fields end
mutable struct LAWaterModel <: AbstractWaterModel @wm_fields end

"Models that use direction variables."
AbstractDirectedModel = Union{CRDWaterModel,LRDWaterModel,QRDWaterModel,CQRDWaterModel}

"Models that don't use direction variables."
AbstractUndirectedModel = Union{NCWaterModel,LAWaterModel}

"Models that include nonlinearities."
AbstractNonlinearModel = Union{NCWaterModel,CRDWaterModel}
