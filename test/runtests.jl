using WaterModels
import InfrastructureModels
import Memento

import Cbc
import Ipopt
import JuMP
import Juniper
import JSON

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

using Test

# Suppress warnings during testing.
Memento.setlevel!(Memento.getlogger(Cbc), "error")
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(WaterModels), "error")

# Default MIP and NLP optimizers.
const cbc = JuMP.with_optimizer(Cbc.Optimizer, logLevel=0)
const ipopt = JuMP.with_optimizer(Ipopt.Optimizer, tol=1.0e-9, acceptable_tol=1.0e-9, max_iter=9999, print_level=0)
const ipopt_ws = JuMP.with_optimizer(Ipopt.Optimizer, tol=1.0e-9, acceptable_tol=1.0e-9, max_iter=9999,
                                     mu_init=1.0e-9, start_with_resto="yes", print_level=0)

include("common.jl")

@testset "WaterModels" begin

    include("base.jl")

    include("data.jl")

    include("multinetwork.jl")

    include("io.jl")

    include("wf.jl")

    include("ne.jl")

    include("owf.jl")

    include("warm_start.jl")

end
