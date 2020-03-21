using WaterModels

import Cbc
import JuMP
import JSON
import InfrastructureModels
import Ipopt
import Juniper
import Memento

# Suppress warnings during testing.
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(Ipopt), "error")
Memento.setlevel!(Memento.getlogger(Juniper), "error")
WaterModels.logger_config!("error")

using Test

# Default MIP and NLP optimizers.
cbc = JuMP.optimizer_with_attributes(Cbc.Optimizer, "logLevel"=>0)
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1.0e-9,
    "acceptable_tol"=>1.0e-9, "max_iter"=>9999, "print_level"=>0)
ipopt_ws = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1.0e-9,
    "acceptable_tol"=>1.0e-9, "max_iter"=>9999, "mu_init"=>1.0e-9,
    "start_with_resto"=>"yes", "print_level"=>0)

include("common.jl")

@testset "WaterModels" begin

    include("base.jl")

    include("data.jl")

    include("io.jl")

    include("wf.jl")

    #include("cwf.jl")

    #include("owf.jl")

    #include("ne.jl")

    #include("warm_start.jl")

end
