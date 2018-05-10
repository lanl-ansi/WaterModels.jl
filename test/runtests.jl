using WaterModels
using Base.Test
using JuMP
using InfrastructureModels
using Ipopt
using Memento
using Pajarito

# Suppress warnings during testing.
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(WaterModels), "error")

# Solver setups.
solver = PajaritoSolver(cont_solver = IpoptSolver())

# Perform the tests.
@testset "WaterModels" begin
    include("data.jl")
    include("feasibility.jl")
end
