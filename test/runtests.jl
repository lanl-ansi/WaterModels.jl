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
const ipopt = JuMP.with_optimizer(Ipopt.Optimizer, tol=1.0e-9, max_iter=9999, print_level=0)
const ipopt_ws = JuMP.with_optimizer(Ipopt.Optimizer, tol=1.0e-9, max_iter=9999, mu_init=1.0e-9, start_with_resto="yes", print_level=0)

@testset "WaterModels" begin

    include("common.jl")

    include("io.jl")

    include("ne.jl")

    include("warm_start.jl")

    include("wf.jl")

end
