using WaterModels
using Base.Test
using JuMP, AmplNLWriter
using InfrastructureModels
using Ipopt
using Memento
using Pajarito

# Suppress warnings during testing.
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(WaterModels), "error")

# Solver setup.
options = ["bonmin.bb_log_level=0","bonmin.fp_log_level=0",
           "bonmin.lp_log_level=0","bonmin.milp_log_level=0",
           "bonmin.nlp_log_level=0","bonmin.oa_log_level=0"]
solver = AmplNLSolver("bonmin", options)

# Perform the tests.
@testset "WaterModels" begin
    include("data.jl")
    include("feasibility.jl")
end
