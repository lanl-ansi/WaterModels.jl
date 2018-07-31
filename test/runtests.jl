using WaterModels
using AmplNLWriter
using Base.Test
using JuMP
using InfrastructureModels
using Memento
using Ipopt

# Suppress warnings during testing.
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(WaterModels), "error")

# Solver setup.
bonmin = AmplNLSolver("bonmin")
ipopt = IpoptSolver()

# Perform the tests.
@testset "WaterModels" begin
    include("data.jl")
    include("feasibility.jl")
end
