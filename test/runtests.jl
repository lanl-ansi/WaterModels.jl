using WaterModels
using GLPK
using GLPKMathProgInterface
using Gurobi
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
glpk = GLPKSolverMIP(msg_lev = GLPK.MSG_ON)
gurobi = GurobiSolver(OutputFlag = 1)
ipopt = IpoptSolver(print_level = 1)
pavito = PavitoSolver(cont_solver = ipopt, mip_solver = glpk)

# Perform the tests.
@testset "WaterModels" begin
    include("data.jl")
    include("hw.jl")
    #include("wf_hw.jl")
    #include("wf_dw.jl")
end
