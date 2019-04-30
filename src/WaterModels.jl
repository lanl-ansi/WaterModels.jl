module WaterModels

import InfrastructureModels
import JSON
import JuMP
import Memento

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

# Create our module level logger (this will get precompiled)
const LOGGER = Memento.getlogger(@__MODULE__)

# Register the module level logger at runtime so that folks can access the logger via `getlogger(WaterModels)`
# NOTE: If this line is not included then the precompiled `WaterModels.LOGGER` won't be registered at runtime.
__init__() = Memento.register(LOGGER)

"Suppresses information and warning messages output by WaterModels. For fine-grained control use the Memento package."
function silence()
    Memento.info(LOGGER, "Suppressing information and warning messages for the rest of this session  Use the Memento package for fine-grained control of logging.")
    Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
    Memento.setlevel!(Memento.getlogger(WaterModels), "error")
end

include("io/common.jl")
include("io/epanet.jl")

include("core/base.jl")
include("core/data.jl")
include("core/function.jl")
include("core/solution.jl")
include("core/types.jl")
include("core/variable.jl")

include("core/constraint.jl")
include("core/constraint_template.jl")
include("core/objective.jl")

include("form/cnlp.jl")
include("form/micp.jl")
include("form/milp_r.jl")
include("form/minlp.jl")
include("form/ncnlp.jl")
include("form/shared.jl")

include("prob/ne.jl")
include("prob/wf.jl")

end
