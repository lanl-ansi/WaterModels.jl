using WaterModels
using Base.Test
using Gurobi
using JuMP, AmplNLWriter
using InfrastructureModels
using Ipopt
using Memento
using Pajarito

# Suppress warnings during testing.
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(WaterModels), "error")

# Solver setup.
# solver = AmplNLSolver("bonmin")
# solver = IpoptSolver()
# solver = PajaritoSolver()
solver = GurobiSolver()

# Perform the tests.
@testset "WaterModels" begin
    #include("data.jl")
    #include("expansion.jl")
    include("feasibility.jl")
end
