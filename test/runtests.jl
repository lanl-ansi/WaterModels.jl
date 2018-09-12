using WaterModels
using AmplNLWriter
using Base.Test
using Cbc
using Gurobi
using JuMP
using InfrastructureModels
using Memento
using Ipopt

# Suppress warnings during testing.
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(WaterModels), "error")

# Solver setup.
cbc = GurobiSolver() #CbcSolver(logLevel = 1)
bonmin = AmplNLSolver("bonmin")
scip = AmplNLSolver("scipampl")
ipopt = IpoptSolver(print_level = 0)

# Perform the tests.
@testset "WaterModels" begin
    #include("data.jl")
    #include("wf_hw.jl")
    include("ne_hw.jl")
    #include("wf_dw.jl")
end
