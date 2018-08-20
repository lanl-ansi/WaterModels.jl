# Define NLP implementations of water distribution models.

export NLPWaterModel, StandardNLPForm

""
@compat abstract type AbstractNLPForm <: AbstractMINLPForm end

""
@compat abstract type StandardNLPForm <: AbstractNLPForm end

"The default NLP model."
const NLPWaterModel = GenericWaterModel{StandardNLPForm}

"Default NLP constructor."
NLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardNLPForm; kwargs...)

"Darcy-Weisbach constraint for flow with unknown direction."
function constraint_dw_unknown_direction{T <: StandardNLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    error(LOGGER, "NLP problem of type wf currently only supports flow variables with known directions.")
end

"Hazen-Williams constraint for flow with unknown direction."
function constraint_hw_unknown_direction{T <: StandardNLPForm}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    error(LOGGER, "NLP problem of type wf currently only supports flow variables with known directions.")
end
