using WaterModels
using Cbc
using GLPKMathProgInterface
using JuMP
using InfrastructureModels
using Ipopt
using Memento
using Pavito
using Test

# Suppress warnings during testing.
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(WaterModels), "error")

# Solver setup.
#cbc = CbcSolver(logLevel = 0, integerTolerance = 1.0e-9,
#                primalTolerance = 1.0e-7, ratioGap = 0.0)
glpk = GLPKSolverMIP(msg_lev = 0, tol_int = 1.0e-9, tol_bnd = 1.0e-7,
                     mip_gap = 0.0, presolve = false)
ipopt = IpoptSolver(print_level = 0)
pavito = PavitoSolver(cont_solver = ipopt, mip_solver = glpk,
                      mip_solver_drives = true)

# Perform the tests.
@testset "WaterModels" begin
    include("data.jl")
    include("wf_hw.jl")
    #include("wf_dw.jl")
    #include("ne_hw.jl")
end
