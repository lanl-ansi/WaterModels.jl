"Undirected nonlinear, nonconvex models."
abstract type AbstractNCModel <: AbstractWaterModel end
mutable struct NCWaterModel <: AbstractNCModel @wm_fields end

"Undirected piecewise linear approximation models."
abstract type AbstractLAModel <: AbstractNCModel end
mutable struct LAWaterModel <: AbstractLAModel @wm_fields end

"Directed nonlinear, nonconvex models."
abstract type AbstractNCDModel <: AbstractNCModel end
mutable struct NCDWaterModel <: AbstractNCDModel @wm_fields end

"Directed nonlinear, quadratic nonconvex models."
abstract type AbstractQRDModel <: AbstractNCDModel end
mutable struct QRDWaterModel <: AbstractQRDModel @wm_fields end

"Directed nonlinear, convex models."
abstract type AbstractCRDModel <: AbstractNCDModel end
mutable struct CRDWaterModel <: AbstractCRDModel @wm_fields end

"Directed nonlinear, convex quadratic models."
abstract type AbstractCQRDModel <: AbstractCRDModel end
mutable struct CQRDWaterModel <: AbstractCQRDModel @wm_fields end

"Directed piecewise-linear relaxation-based models."
abstract type AbstractPWLRDModel <: AbstractCQRDModel end
mutable struct PWLRDWaterModel <: AbstractPWLRDModel @wm_fields end

"Directed linear relaxation-based models."
abstract type AbstractLRDModel <: AbstractCQRDModel end
mutable struct LRDWaterModel <: AbstractLRDModel @wm_fields end

"Models that include nonquadratic nonlinearities."
AbstractNonlinearModel = Union{NCWaterModel, NCDWaterModel, CRDWaterModel}