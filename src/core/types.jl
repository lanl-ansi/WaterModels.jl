@enum FLOW_DIRECTION POSITIVE=1 NEGATIVE=-1 UNKNOWN=0

export CNLPWaterModel, StandardCNLPForm,
    NCNLPWaterModel, StandardNCNLPForm,
    MILPRWaterModel, StandardMILPRForm,
    MICPWaterModel, StandardMICPForm


"AbstractCNLPForm is derived from AbstractWaterFormulation"
abstract type AbstractCNLPForm <: AbstractWaterFormulation end

"StandardCNLPForm is derived from AbstractCNLPForm"
abstract type StandardCNLPForm <: AbstractCNLPForm end

"The CNLP model relaxes constraints into the objective."
const CNLPWaterModel = GenericWaterModel{StandardCNLPForm}

"Default CNLP constructor"
CNLPWaterModel(data::Dict{String, Any}; kwargs...) = GenericWaterModel(data, StandardCNLPForm; kwargs...)


abstract type AbstractNCNLPForm <: AbstractWaterFormulation end
abstract type StandardNCNLPForm <: AbstractNCNLPForm end

"The default NCNLP (non-convex nonlinear programming) model retains the exact head loss physics."
const NCNLPWaterModel = GenericWaterModel{StandardNCNLPForm}

"Default NCNLP constructor"
NCNLPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardNCNLPForm; kwargs...)



abstract type AbstractMILPRForm <: AbstractWaterFormulation end
abstract type StandardMILPRForm <: AbstractMILPRForm end

"The default MILPR (mixed-integer linear, relaxed) model is a linear outer-approximation of the MICP model."
const MILPRWaterModel = GenericWaterModel{StandardMILPRForm}

"Default MILPR constructor"
MILPRWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMILPRForm; kwargs...)



abstract type AbstractMICPForm <: AbstractWaterFormulation end
abstract type StandardMICPForm <: AbstractMICPForm end

"The default MICP (mixed-integer convex program) model is a relaxation of the non-convex MINLP model."
const MICPWaterModel = GenericWaterModel{StandardMICPForm}

"Default MICP constructor."
MICPWaterModel(data::Dict{String,Any}; kwargs...) = GenericWaterModel(data, StandardMICPForm; kwargs...)
