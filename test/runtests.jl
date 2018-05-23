using WaterModels
using Base.Test
using JuMP, AmplNLWriter
using InfrastructureModels
using Ipopt
using Memento

# Suppress warnings during testing.
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(WaterModels), "error")

# Solver setup.
#solver = AmplNLSolver("couenne")
solver = IpoptSolver()

# Perform the tests.
@testset "WaterModels" begin
    include("data.jl")
    include("feasibility.jl")
end