using WaterModels
using AmplNLWriter
using Base.Test
using Gurobi
using JuMP
using InfrastructureModels
using Memento
using Pajarito

# Suppress warnings during testing.
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(WaterModels), "error")

# Solver setup.
# solver = AmplNLSolver("bonmin")
# solver = GurobiSolver(FeasibilityTol = 1.0e-9, OptimalityTol = 1.0e-9, MIPGap = 0.0)
# solver = PajaritoSolver(mip_solver = GurobiSolver())
solver = GurobiSolver()

# Perform the tests.
@testset "WaterModels" begin
    include("data.jl")
    #include("expansion.jl")
    include("feasibility.jl")
end
