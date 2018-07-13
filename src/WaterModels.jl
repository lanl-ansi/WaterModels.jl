isdefined(Base, :__precompile__) && __precompile__()

module WaterModels

using JSON
using InfrastructureModels
using MathProgBase
using JuMP
using Compat
using Memento

import Compat: @__MODULE__

# Create and initialize the module level logger.
const LOGGER = getlogger(@__MODULE__)
__init__() = Memento.register(LOGGER)

include("io/common.jl")
include("io/epanet.jl")

include("core/base.jl")
include("core/constraint.jl")
include("core/data.jl")
include("core/objective.jl")
include("core/solution.jl")
include("core/variable.jl")

include("form/minlp.jl")

include("prob/expansion.jl")
include("prob/feasibility.jl")

end
