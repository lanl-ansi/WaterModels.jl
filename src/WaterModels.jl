isdefined(Base, :__precompile__) && __precompile__()

module WaterModels

using JSON
using InfrastructureModels
using MathProgBase
using JuMP
using Compat
using Memento

import Compat: @__MODULE__

# Create the module level logger. (This will get precompiled.)
const LOGGER = getlogger(@__MODULE__)
setlevel!(LOGGER, "info")

# Register the module level logger at runtime so the logger can be accessed via `getlogger(PowerModels)`
# NOTE: If this line is not included, then the precompiled `WaterModels.LOGGER` won't be registered at runtime.
__init__() = Memento.register(LOGGER)

include("io/epanet.jl")
include("io/json.jl")
include("io/common.jl")

include("core/data.jl")
include("core/base.jl")
include("core/variable.jl")
include("core/constraint.jl")
include("core/objective.jl")
include("core/solution.jl")

include("form/minlp.jl")
include("form/misocp.jl")

include("prob/wf.jl")

end
