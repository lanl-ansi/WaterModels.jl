using WaterModels
import InfrastructureModels
import Memento

# Suppress warnings during testing.
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(WaterModels), "error")

import Cbc
import Ipopt
import Juniper
import JuMP
import JSON

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

using Test

# Default MIP, NLP, and MINLP optimizers.
const cbc = JuMP.with_optimizer(Cbc.Optimizer, logLevel=0)
const ipopt = JuMP.with_optimizer(Ipopt.Optimizer, tol=1.0e-9, print_level=0)
const juniper = JuMP.with_optimizer(Juniper.Optimizer, nl_solver=ipopt, log_levels=[])

@testset "WaterModels" begin

    include("io.jl")

    include("ne.jl")

    include("wf.jl")

end
