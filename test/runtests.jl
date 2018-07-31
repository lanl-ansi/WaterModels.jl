using WaterModels
using Base.Test
using Cbc
using JuMP
using InfrastructureModels
using Memento
using Ipopt
using Pavito

# Suppress warnings during testing.
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(WaterModels), "error")

# Solver setup.
cbc = CbcSolver(logLevel = 0)
ipopt = IpoptSolver(print_level = 0)
pavito = PavitoSolver(mip_solver = cbc, cont_solver = ipopt,
                      mip_solver_drives = false, log_level = 0)

# Perform the tests.
@testset "WaterModels" begin
    #include("data.jl")
    include("feasibility.jl")
end
