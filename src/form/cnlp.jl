# Define CNLP implementations of water distribution models.

export CNLPWaterModel, StandardCNLPForm

import MathProgBase

"AbstractCNLPForm is derived from AbstractWaterFormulation"
abstract type AbstractCNLPForm <: AbstractWaterFormulation end

"StandardCNLPForm is derived from AbstractCNLPForm"
abstract type StandardCNLPForm <: AbstractCNLPForm end

"The CNLP model is a convex formulation used for obtaining flow rates for a fixed network design."
const CNLPWaterModel = GenericWaterModel{StandardCNLPForm}

"CNLP constructor."
CNLPWaterModel(data::Dict{String, Any}; kwargs...) = GenericWaterModel(data, StandardCNLPForm; kwargs...)

function variable_directed_flow(wm::GenericWaterModel{T}, n::Int = wm.cnw) where T <: StandardCNLPForm
    # Get indices for all network arcs.
    arcs = collect(ids(wm, n, :connection))

    # Initialize directed flow variables. The variables qp correspond to flow
    # from i to j, and the variables qn correspond to flow from j to i.
    wm.var[:nw][n][:qp] = Dict{Int, Array{Variable, 1}}()
    wm.var[:nw][n][:qn] = Dict{Int, Array{Variable, 1}}()

    for (a, connection) in wm.ref[:nw][n][:connection]
        # Get the number of possible resistances for this arc.
        num_resistances = length(wm.ref[:nw][n][:resistance][a])

        # Initialize variables associated with flow from i to j.
        wm.var[:nw][n][:qp][a] = @variable(wm.model, [r in 1:num_resistances],
                                           lowerbound = 0.0, start = 0.0,
                                           category = :Cont,
                                           basename = "qp_$(n)_$(a)")

        # Initialize variables associated with flow from j to i.
        wm.var[:nw][n][:qn][a] = @variable(wm.model, [r in 1:num_resistances],
                                           lowerbound = 0.0, start = 0.0,
                                           category = :Cont,
                                           basename = "qn_$(n)_$(a)")
    end
end
