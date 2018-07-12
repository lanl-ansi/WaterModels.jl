# Define the MINLP implementations of water distribution models.

export MINLPWaterModel, StandardMINLPForm

abstract type AbstractMINLPForm <: AbstractWaterFormulation end

abstract type StandardMINLPForm <: AbstractMINLPForm end

const MINLPWaterModel = GenericWaterModel{StandardMINLPForm}

MINLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMINLPForm; kwargs...)
