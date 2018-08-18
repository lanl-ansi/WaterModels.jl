# Define NLP implementations of water distribution models.

export NLPWaterModel, StandardNLPForm

""
@compat abstract type AbstractNLPForm <: AbstractMINLPForm end

""
@compat abstract type StandardNLPForm <: AbstractNLPForm end

const NLPWaterModel = GenericWaterModel{StandardNLPForm}

"Default NLP constructor."
NLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardNLPForm; kwargs...)
