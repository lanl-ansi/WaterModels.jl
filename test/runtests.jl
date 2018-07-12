using WaterModels
using AmplNLWriter
using Base.Test
using Gurobi
using JuMP
using InfrastructureModels
using Memento

# Suppress warnings during testing.
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(WaterModels), "error")

# Solver setup.
solver = GurobiSolver()
# solver = AmplNLSolver("bonmin")

# Perform the tests.
@testset "WaterModels" begin
    include("data.jl")
    include("feasibility.jl")
end
