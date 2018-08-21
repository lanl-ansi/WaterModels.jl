# Define NLP implementations of water distribution models.

export NLPWaterModel, StandardNLPForm

""
@compat abstract type AbstractNLPForm <: AbstractMINLPForm end

""
@compat abstract type StandardNLPForm <: StandardMINLPForm end

"The default NLP model."
const NLPWaterModel = GenericWaterModel{StandardNLPForm}

"Default NLP constructor."
NLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardNLPForm; kwargs...)

"Hazen-Williams constraint for flow with unknown direction."
function constraint_hw_unknown_direction{T <: StandardNLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    error(LOGGER, "NLP problem of type wf_hw currently only supports flow variables with known directions.")
end
