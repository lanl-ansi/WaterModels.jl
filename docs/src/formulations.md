# Network Formulations
```@meta
CurrentModule = WaterModels
```

All methods for constructing WaterModels are ultimately derived from the following type:
```@docs
AbstractWaterModel
```

## Type Hierarchy
We begin with the top of the modeling hierarchy, where we can distinguish among approximations and relaxations for the the physics of water flow.
There are currently six formulations supported in WaterModels.
Specifically, there are two exact nonconvex formulations, one approximation, and three relaxations:
```julia
# Exact formulations.
abstract type AbstractNCModel <: AbstractWaterModel end
abstract type AbstractNCDModel <: AbstractNCModel end

# Approximation-based formulations.
abstract type AbstractLAModel <: AbstractNCModel end

# Relaxation-based formulations.
abstract type AbstractCRDModel <: AbstractNCDModel end
abstract type AbstractLRDModel <: AbstractCRDModel end
abstract type AbstractPWLRDModel <: AbstractCRDModel end
```

## Water Models
Each of the abstract modeling types are used to derive a more specific WaterModel type:
```julia
mutable struct NCWaterModel <: AbstractNCModel @wm_fields end
mutable struct NCDWaterModel <: AbstractNCDModel @wm_fields end
mutable struct LAWaterModel <: AbstractLAModel @wm_fields end
mutable struct CRDWaterModel <: AbstractCRDModel @wm_fields end
mutable struct LRDWaterModel <: AbstractLRDModel @wm_fields end
mutable struct PWLRDWaterModel <: AbstractPWLRDModel @wm_fields end
```

## User-Defined Abstractions
User-defined abstractions can begin from a root abstract model like `AbstractWaterModel`, e.g.,
```julia
abstract type AbstractFooModel <: AbstractWaterModel end
mutable struct FooWaterModel <: AbstractFooModel @wm_fields end
```
They can also be derived from existing modeling types, e.g.,
```julia
abstract type AbstractSpecialWaterModel <: AbstractPWLRDModel end
mutable struct SpecialWaterModel <: AbstractSpecialWaterModel @wm_fields end
```
As an example, the above might be used to define some specialized piecewise-linear, direction-based water flow relaxation.
This definition of new types (`AbstractSpecialWaterModel` and `SpecialWaterModel`) would allow the developer to reuse common variables, constraints, and objectives already defined over `AbstractPWLRDModel` while also potentially overriding other variables, constraints, and objectives.
This overriding is accomplished by implementing new but equivalently-named model-building functions that are defined over `AbstractSpecialWaterModel` instead of `AbstractPWLRDModel` (or its ancestors, i.e., `AbstractCRDModel`, `AbstractNCDModel`, and `AbstractWaterModel`).

## Supported Formulations
All formulation names refer to how the underlying physics of a water network are modeled.
For example, the `LRD` model uses a linear relaxation (LR), direction-based (D) model of water network physics.

| Formulation      | Description           |
| ---------------- | --------------------- |
| NC               | Physics modeled using nonlinear equations. |
| NCD              | Physics modeled using nonlinear equations. Flow direction modeled using binary variables. |
| LA               | Physics modeled using piecewise-linear approximations. |
| CRD              | Physics modeled using convex relaxations. Flow direction modeled using binary variables. |
| LRD              | Physics modeled using linear relaxations. Flow direction modeled using binary variables. |
| PWLRD            | Physics modeled using piecewise relaxations. Flow direction modeled using binary variables. |
