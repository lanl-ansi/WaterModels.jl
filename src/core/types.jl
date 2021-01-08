"Undirected nonlinear, nonconvex models."
abstract type AbstractNCModel <: AbstractWaterModel end
mutable struct NCWaterModel <: AbstractNCModel @wm_fields end

"Undirected piecewise linear approximation models."
abstract type AbstractLAModel <: AbstractNCModel end
mutable struct LAWaterModel <: AbstractLAModel @wm_fields end

"Directed nonlinear, nonconvex models."
abstract type AbstractNCDModel <: AbstractNCModel end
mutable struct NCDWaterModel <: AbstractNCDModel @wm_fields end

"Directed nonlinear, convex models."
abstract type AbstractCRDModel <: AbstractNCDModel end
mutable struct CRDWaterModel <: AbstractCRDModel @wm_fields end

"Directed piecewise-linear relaxation-based models."
abstract type AbstractPWLRDModel <: AbstractCRDModel end
mutable struct PWLRDWaterModel <: AbstractPWLRDModel @wm_fields end

"Directed linear relaxation-based models."
abstract type AbstractLRDModel <: AbstractCRDModel end
mutable struct LRDWaterModel <: AbstractLRDModel @wm_fields end

"Models that include nonquadratic nonlinearities."
AbstractNonlinearModel = Union{NCWaterModel, NCDWaterModel, CRDWaterModel}