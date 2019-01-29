module WaterModels

using JSON
using InfrastructureModels
using MathProgBase
using JuMP
using Memento

# Create and initialize the module level logger.
const LOGGER = getlogger(@__MODULE__)
__init__() = Memento.register(LOGGER)

include("io/common.jl")
include("io/epanet.jl")

include("core/types.jl")

include("core/base.jl")
include("core/constraint.jl")
include("core/constraint_template.jl")
include("core/data.jl")
include("core/function.jl")
include("core/objective.jl")
include("core/solution.jl")
include("core/variable.jl")

include("form/cvxnlp.jl")
include("form/minlp.jl")
include("form/milp_r.jl")
include("form/shared.jl")
#include("form/micp.jl")
#include("form/milp.jl")
#include("form/milp_r.jl")
#include("form/minlp_b.jl")
#include("form/nlp.jl")
#include("form/equality.jl")
#include("form/exact.jl")
#include("form/relaxed.jl")

include("prob/cvx.jl")
include("prob/ne.jl")
include("prob/wf.jl")

include("util/bound_tightening.jl")
include("util/get_cvx_solution.jl")
include("util/compute_objective.jl")
include("util/repair_solution.jl")
include("util/find_initial_solution.jl")
include("util/set_initial_solution.jl")
include("util/eliminate_variables.jl")
include("util/user_cuts.jl")
include("util/bounding_cuts.jl")
include("util/lazy_cuts.jl")
include("util/heuristic_cuts.jl")
include("util/solve_global.jl")

end
