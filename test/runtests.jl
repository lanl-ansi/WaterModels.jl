using WaterModels
import InfrastructureModels
import Memento

# Suppress warnings during testing.
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(WaterModels), "error")

import Cbc
import Ipopt
import JuMP
import JSON

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

using Test

# Default MIP and NLP optimizers.
const cbc = JuMP.with_optimizer(Cbc.Optimizer, logLevel=0)
const ipopt = JuMP.with_optimizer(Ipopt.Optimizer, tol=1.0e-9, print_level=0)

@testset "WaterModels" begin

    include("base.jl")

    include("io.jl")

    include("ne.jl")

    include("wf.jl")

end
