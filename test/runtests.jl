using WaterModels
using InfrastructureModels
using Memento

# Suppress warnings during testing.
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(WaterModels), "error")

using JuMP, AmplNLWriter
using Base.Test

# Solver setups.
bonmin = AmplNLSolver("bonmin", ["bonmin.nlp_log_level=0"])

# Perform the tests.
@testset "WaterModels" begin

include("data.jl")
include("feasibility.jl")

end
