using WaterModels
using GLPKMathProgInterface
using InfrastructureModels
using Ipopt
using JuMP
using Memento
using Pavito
using Test

# Suppress warnings during testing.
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(WaterModels), "error")

# Solver setup.
glpk = GLPKSolverMIP(presolve = false)
ipopt = IpoptSolver(print_level = 1, tol = 1.0e-9)
pavito = PavitoSolver(cont_solver = ipopt, mip_solver = glpk)

# Perform the tests.
@testset "WaterModels" begin
    include("data.jl")
    include("hw.jl")
end
