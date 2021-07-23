using WaterModels: _calc_flow_per_unit_transform, _calc_head_per_unit_transform
using WaterModels

import Cbc
import JuMP
import JSON
import Ipopt
import Juniper
import Logging
import Memento

const _IM = WaterModels._IM
const _MOI = WaterModels._MOI

# Suppress warnings during testing.
Memento.setlevel!(Memento.getlogger(_IM), "error")
Memento.setlevel!(Memento.getlogger(Ipopt), "error")
Memento.setlevel!(Memento.getlogger(Juniper), "error")
WaterModels.logger_config!("error")
Logging.disable_logging(Logging.Info)

using Test

# Default MIP optimizer
cbc = JuMP.optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)

# Default NLP optimizer.
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1.0e-6, "print_level" => 0, "sb" => "yes")

include("common.jl")

@testset "WaterModels" begin

    include("epanet.jl")

    include("base.jl")

    include("node.jl")

    include("demand.jl")

    include("reservoir.jl")

    include("tank.jl")

    include("pipe.jl")

    include("pump.jl")

    include("data.jl")

    include("function.jl")

    include("io.jl")

    include("wf.jl")

    # include("owf.jl")

    # include("des.jl")

    # include("relax.jl")

    # include("variable_index.jl")

    # include("pairwise_cuts.jl")

    # include("obbt.jl")

    # include("warm_start.jl")

end
