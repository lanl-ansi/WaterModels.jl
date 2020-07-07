module WaterModels
    import InfrastructureModels
    import InfrastructureModels: ids, ref, var, con, sol, nw_ids, nws, optimize_model!, @im_fields
    const _IM = InfrastructureModels

    import JSON
    import JuMP
    import LsqFit
    import Memento

    import MathOptInterface
    const _MOI = MathOptInterface

    # Create our module-level logger (this will get precompiled).
    const _LOGGER = Memento.getlogger(@__MODULE__)

    # Register the module-level logger at runtime so users can access the logger
    # via `getlogger(WaterModels)` NOTE: If this line is not included, then the
    # precompiled `WaterModels._LOGGER` will not be registered at runtime.
    __init__() = Memento.register(_LOGGER)

    "Suppresses information and warning messages output by WaterModels. For
    more fine-grained control, use the Memento package."
    function silence()
        Memento.info(_LOGGER, "Suppressing information and warning messages for "
            * "the rest of this session. Use the Memento package for more "
            * "fine-grained control of logging.")
        Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
        Memento.setlevel!(Memento.getlogger(WaterModels), "error")
    end

    "Allows the user to set the logging level without the need to add Memento."
    function logger_config!(level)
        Memento.config!(Memento.getlogger("WaterModels"), level)
    end

    const _wm_global_keys = Set(["time_series", "per_unit", "option"])

    include("io/common.jl")
    include("io/epanet.jl")

    include("core/base.jl")
    include("core/data.jl")
    include("core/ref.jl")
    include("core/types.jl")
    include("core/function.jl")
    include("core/variable.jl")

    include("core/constraint.jl")
    include("core/constraint_template.jl")
    include("core/objective.jl")

    include("form/directed_flow.jl")
    include("form/undirected_flow.jl")
    include("form/micp.jl")
    include("form/milp.jl")
    include("form/milpr.jl")
    include("form/nlp.jl")

    include("prob/wf.jl")
    include("prob/owf.jl")
    include("prob/des.jl")

    include("core/export.jl")
end
