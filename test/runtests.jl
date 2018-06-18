using WaterModels
using Base.Test
using DataFrames
using JuMP, AmplNLWriter
using InfrastructureModels
using Memento

# Suppress warnings during testing.
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(WaterModels), "error")

# Solver setup.
solver = AmplNLSolver("bonmin")

# Perform the tests.
@testset "WaterModels" begin
    include("data.jl")
    include("feasibility.jl")
end
