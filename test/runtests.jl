using WaterModels

import InfrastructureModels
const _IM = InfrastructureModels

import Cbc
import JuMP
import JSON
import Ipopt
import Juniper
import Memento
import MathOptInterface
const _MOI = MathOptInterface

# Suppress warnings during testing.
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(Ipopt), "error")
Memento.setlevel!(Memento.getlogger(Juniper), "error")
WaterModels.logger_config!("error")

using Test

# Default MIP and NLP optimizers.
cbc = JuMP.optimizer_with_attributes(Cbc.Optimizer, "logLevel"=>0)

ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1.0e-9,
    "acceptable_tol"=>1.0e-9, "max_iter"=>9999, "print_level"=>0, "sb"=>"yes")

ipopt_ws = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1.0e-9,
    "acceptable_tol"=>1.0e-9, "max_iter"=>9999, "mu_init"=>1.0e-9,
    "start_with_resto"=>"yes", "print_level"=>0, "sb"=>"yes")

include("common.jl")

@testset "WaterModels" begin

    include("base.jl")

    include("data.jl")

    include("io.jl")

    include("wf.jl")

    include("owf.jl")

    include("owf_agm.jl")

    include("des.jl")

    include("warm_start.jl")

end
