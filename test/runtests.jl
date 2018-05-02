using WaterModels
using InfrastructureModels
using Memento

# Suppress warnings during testing.
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(WaterModels), "error")

using Ipopt
using Base.Test

# Solver setups.
ipopt_solver = IpoptSolver(tol = 1.0e-6, print_level = 0)

# Perform the tests.
@testset "WaterModels" begin

include("data.jl")

end
