isdefined(Base, :__precompile__) && __precompile__()

module WaterModels

using JSON
using MathProgBase
using JuMP
using Compat

#include("io/epanet.jl")
include("io/json.jl")
include("io/common.jl")



include("core/base.jl")
include("core/data.jl")
include("core/variable.jl")
include("core/constraint.jl")
# include("core/objective.jl")
include("core/solution.jl")

include("form/minlp.jl")
include("form/misocp.jl")

include("prob/wf.jl")
# include("prob/ne.jl")
# include("prob/ls.jl")
# include("prob/nels.jl")
# include("prob/nelsfd.jl")

end
