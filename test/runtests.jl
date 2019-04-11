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

using Test

# default setup for optimizers
ipopt = JuMP.with_optimizer(Ipopt.Optimizer, tol=1.0e-9, print_level=5)
cbc = JuMP.with_optimizer(Cbc.Optimizer, logLevel=1)

@testset "WaterModels" begin

    #include("io.jl")

    #include("ne.jl")

    include("wf.jl")

end
