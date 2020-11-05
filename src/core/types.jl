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
