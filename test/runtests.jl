using WaterModels
using CPLEX
using GLPK
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
#glpk = GLPKSolverMIP(msg_lev = 0, tol_int = 1.0e-9, tol_bnd = 1.0e-7, mip_gap = 0.0, presolve = false)
glpk = CplexSolver(CPX_PARAM_EPOPT = 1.0e-6)
#glpk = GLPKSolverMIP(msg_lev = GLPK.MSG_ON, tol_int = 1.0e-9, tol_bnd = 1.0e-7, mip_gap = 0.0, presolve = false)
ipopt = IpoptSolver(print_level = 0)
pavito = PavitoSolver(cont_solver = ipopt, mip_solver = glpk)

# Perform the tests.
@testset "WaterModels" begin
    include("data.jl")
    include("hw.jl")
    #include("wf_hw.jl")
    #include("wf_dw.jl")
end
