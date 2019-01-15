# WaterModels.jl Documentation
```@meta
CurrentModule = WaterModels
```

## Overview
WaterModels.jl is a Julia package for steady state water network optimization.
It is designed to enable computational evaluation of historical and emerging water network formulations and algorithms using a common platform.
The software is engineered to decouple problem specifications (e.g., feasibility, network expansion) from water network optimization formulations (e.g., mixed-integer linear, mixed-integer nonlinear).
This decoupling enables the definition of a wide variety of water network optimization formulations and their comparison on common problem specifications.

## Installation
The latest stable release of WaterModels can be installed using the Julia package manager with

```julia
using Pkg
Pkg.add("WaterModels")
```

For the current development version, "check out" this package with
```julia
using Pkg
Pkg.checkout("WaterModels")
```

You can test that the package works by running
```julia
using Pkg
Pkg.test("WaterModels")
```
