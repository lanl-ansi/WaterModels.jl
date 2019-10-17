# WaterModels.jl Documentation

```@meta
CurrentModule = WaterModels
```

## Overview
WaterModels.jl is a Julia/JuMP package for steady state water network optimization.
It is designed to enable computational evaluation of historical and emerging water network formulations and algorithms using a common platform.
The code is engineered to decouple [Problem Specifications](@ref) (e.g., water flow, optimal water flow, network expansion) from [Network Formulations](@ref) (e.g., mixed-integer linear, mixed-integer nonlinear).
This decoupling enables the definition of a wide variety of water network optimization formulations and their comparison on common problem specifications.

## Installation
The latest stable release of WaterModels can be installed using the Julia package manager with
```julia
] add WaterModels
```

For the current development version, install the package using
```julia
] add WaterModels#master
```

Test that the package works by executing
```julia
] test WaterModels
```
