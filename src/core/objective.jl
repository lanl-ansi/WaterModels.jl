function objective_minimize_gamma(wm::GenericWaterModel)
    arcs_from = collect(ids(wm, :pipes))
    return @objective(wm.model, Min, sum(sum(wm.var[:nw][n][:gamma][a] for a in arcs_from) for n in nws(wm)))
end
