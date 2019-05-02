export
    GenericWaterModel,
    optimize!,
    run_generic_model, build_generic_model, solve_generic_model,
    ismultinetwork, nw_ids, nws,
    ids, ref, var, con, ext

""
abstract type AbstractWaterFormulation end

"""
```
type GenericWaterModel{T<:AbstractWaterFormulation}
    model::JuMP.Model
    data::Dict{String,<:Any}
    setting::Dict{String,<:Any}
    solution::Dict{String,<:Any}
    ref::Dict{Symbol,<:Any} # reference data
    var::Dict{Symbol,<:Any} # JuMP variables
    con::Dict{Symbol,<:Any} # JuMP constraint references
    cnw::Int                # current network index value
    ext::Dict{Symbol,<:Any} # user extentions
end
```
where

* `data` is the original data, usually from reading in a `.inp` file,
* `setting` usually looks something like `Dict("output" => Dict("flows" => true))`, and
* `ref` is a place to store commonly-used precomputed data from the data dictionary,
  primarily for converting datatypes, filtering deactivated components, and storing
  system-wide values that need to be computed globally. See `build_ref(data)` for further details.

Methods on `GenericWaterModel` for defining variables and adding constraints should

* work with the `ref` dict, rather than the original `data` dict,
* add them to `model::JuMP.Model`, and
* follow the conventions for variable and constraint names.
"""
mutable struct GenericWaterModel{T<:AbstractWaterFormulation}
    model::JuMP.Model

    data::Dict{String,<:Any}
    setting::Dict{String,<:Any}
    solution::Dict{String,<:Any}

    ref::Dict{Symbol,<:Any}
    var::Dict{Symbol,<:Any}
    con::Dict{Symbol,<:Any}
    cnw::Int

    # Extensions should define a type to hold information particular to
    # their functionality, and store an instance of the type in this
    # dictionary keyed on an extension-specific symbol.
    ext::Dict{Symbol,<:Any}
end

function GenericWaterModel(data::Dict{String,<:Any}, T::DataType; ext = Dict{Symbol,Any}(), setting = Dict{String,Any}(), jump_model::JuMP.Model=JuMP.Model(), kwargs...)
    ref = build_ref(data)
    var = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    con = Dict{Symbol,Any}(:nw => Dict{Int,Any}())

    for nw_id in keys(ref[:nw])
        var[:nw][nw_id] = Dict{Symbol, Any}()
        con[:nw][nw_id] = Dict{Symbol, Any}()
    end

    cnw = minimum([k for k in keys(var[:nw])])

    wm = GenericWaterModel{T}(
        jump_model,
        data,
        setting,
        Dict{String,Any}(), # solution
        ref,
        var,
        con,
        cnw,
        ext
       )

    return wm
end

### Helper functions for working with multinetworks.
""
ismultinetwork(wm::GenericWaterModel) = (length(wm.ref[:nw]) > 1)

""
nw_ids(wm::GenericWaterModel) = keys(wm.ref[:nw])

""
nws(wm::GenericWaterModel) = wm.ref[:nw]

""
ids(wm::GenericWaterModel, nw::Int, key::Symbol) = keys(wm.ref[:nw][nw][key])
ids(wm::GenericWaterModel, key::Symbol; nw::Int=wm.cnw) = keys(wm.ref[:nw][nw][key])

""
ref(wm::GenericWaterModel, nw::Int) = wm.ref[:nw][nw]
ref(wm::GenericWaterModel, nw::Int, key::Symbol) = wm.ref[:nw][nw][key]
ref(wm::GenericWaterModel, nw::Int, key::Symbol, idx) = wm.ref[:nw][nw][key][idx]
ref(wm::GenericWaterModel, nw::Int, key::Symbol, idx, param::String) = wm.ref[:nw][nw][key][idx][param]

ref(wm::GenericWaterModel; nw::Int=wm.cnw) = wm.ref[:nw][nw]
ref(wm::GenericWaterModel, key::Symbol; nw::Int=wm.cnw) = wm.ref[:nw][nw][key]
ref(wm::GenericWaterModel, key::Symbol, idx; nw::Int=wm.cnw) = wm.ref[:nw][nw][key][idx]

""
var(wm::GenericWaterModel, nw::Int) = wm.var[:nw][nw]
var(wm::GenericWaterModel, nw::Int, key::Symbol) = wm.var[:nw][nw][key]
var(wm::GenericWaterModel, nw::Int, key::Symbol, idx) = wm.var[:nw][nw][key][idx]

""
con(wm::GenericWaterModel, nw::Int) = wm.con[:nw][nw]
con(wm::GenericWaterModel, nw::Int, key::Symbol) = wm.con[:nw][nw][key]
con(wm::GenericWaterModel, nw::Int, key::Symbol, idx) = wm.con[:nw][nw][key][idx]

""
function optimize!(wm::GenericWaterModel, optimizer::JuMP.OptimizerFactory)
    if wm.model.moi_backend.state == MOIU.NO_OPTIMIZER
        _, solve_time, solve_bytes_alloc, sec_in_gc = @timed JuMP.optimize!(wm.model, optimizer)
    else
        Memento.warn(LOGGER, "Model already contains optimizer factory, cannot use optimizer specified in `solve_generic_model`")
        _, solve_time, solve_bytes_alloc, sec_in_gc = @timed JuMP.optimize!(wm.model)
    end

    try
        solve_time = MOI.get(wm.model, MOI.SolveTime())
    catch
        Memento.warn(LOGGER, "The given optimizer does not provide the SolveTime() attribute, falling back on @timed. This is not a rigorous timing value.");
    end

    return JuMP.termination_status(wm.model), JuMP.primal_status(wm.model), JuMP.dual_status(wm.model), solve_time
end

""
function run_generic_model(file::String, model_constructor, optimizer, post_method; relaxed::Bool=false, kwargs...)
    data = WaterModels.parse_file(file)
    return run_generic_model(data, model_constructor, optimizer, post_method, relaxed=relaxed; kwargs...)
end

""
function run_generic_model(data::Dict{String,<:Any}, model_constructor, optimizer, post_method; relaxed::Bool=false, solution_builder=get_solution, kwargs...)
    wm = build_generic_model(data, model_constructor, post_method; kwargs...)
    #wm, time, bytes_alloc, sec_in_gc = @timed build_generic_model(data, model_constructor, post_method; kwargs...)
    #println("model build time: $(time)")

    solution = solve_generic_model(wm, optimizer, relaxed=relaxed; solution_builder = solution_builder)
    #solution, time, bytes_alloc, sec_in_gc = @timed solve_generic_model(wm, optimizer; solution_builder = solution_builder)
    #println("solution time: $(time)")

    return solution
end

""
function build_generic_model(file::String, model_constructor, post_method; kwargs...)
    data = WaterModels.parse_file(file)
    return build_generic_model(data, model_constructor, post_method; kwargs...)
end

""
function build_generic_model(data::Dict{String,<:Any}, model_constructor, post_method; multinetwork=false, kwargs...)
    # NOTE, this model constructor will build the ref dict using the latest info from the data
    wm = model_constructor(data; kwargs...)

    if !multinetwork && ismultinetwork(wm)
        Memento.error(LOGGER, "Attempted to build a single-network model with multi-network data.")
    end

    post_method(wm)

    return wm
end

""
function parse_status(termination_status::MOI.TerminationStatusCode, primal_status::MOI.ResultStatusCode, dual_status::MOI.ResultStatusCode)
    if termination_status == MOI.OPTIMAL
        return :Optimal
    elseif termination_status == MOI.LOCALLY_SOLVED
        return :LocalOptimal
    elseif termination_status == MOI.INFEASIBLE
        return :Infeasible
    elseif termination_status == MOI.LOCALLY_INFEASIBLE
        return :LocalInfeasible
    else
        return :Error
    end
end

""
function solve_generic_model(wm::GenericWaterModel, optimizer::JuMP.OptimizerFactory; relaxed::Bool=false, solution_builder = get_solution)
    if relaxed
        relax_integrality!(wm)
    end

    termination_status, primal_status, dual_status, solve_time = optimize!(wm, optimizer)
    status = parse_status(termination_status, primal_status, dual_status)

    solution = build_solution(wm, status, solve_time; solution_builder = solution_builder)
    #solution, time, bytes_alloc, sec_in_gc = @timed build_solution(wm, status, solve_time; solution_builder = solution_builder)
    #println("build_solution time: $(time)")

    return solution
end

"""
Returns a dict that stores commonly-used, precomputed data from of the data
dictionary, primarily for converting data types, filtering out deactivated
components, and storing system-wide values that need to be computed globally.

Some of the common keys include:

* `:pipes` -- the set of pipes in the network,
* `:pumps` -- the set of pumps in the network,
* `:valves` -- the set of valves in the network,
* `:links` -- the set of all links in the network,
* `:junctions` -- the set of junctions in the network,
* `:reservoirs` -- the set of reservoirs in the network,
* `:tanks` -- the set of tanks in the network,
* `:emitters` -- the set of emitters in the network,
* `:nodes` -- the set of all nodes in the network
"""
function build_ref(data::Dict{String,<:Any})
    refs = Dict{Symbol,Any}()

    nws = refs[:nw] = Dict{Int,Any}()

    if InfrastructureModels.ismultinetwork(data)
        nws_data = data["nw"]
    else
        nws_data = Dict("0" => data)
    end

    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        ref = nws[nw_id] = Dict{Symbol,Any}()

        for (key, item) in nw_data
            if isa(item, Dict{String,Any})
                try
                    item_lookup = Dict{Int,Any}([(parse(Int, k), v) for (k,v) in item])
                    ref[Symbol(key)] = item_lookup
                catch
                    ref[Symbol(key)] = item
                end
            else
                ref[Symbol(key)] = item
            end
        end

        ref[:links] = merge(ref[:pipes], ref[:valves], ref[:pumps])
        ref[:links_ne] = filter(is_ne_link, ref[:links])
        ref[:links_known_direction] = filter(has_known_flow_direction, ref[:links])
        ref[:links_unknown_direction] = filter(!has_known_flow_direction, ref[:links])
        ref[:nodes] = merge(ref[:junctions], ref[:reservoirs], ref[:emitters])

        # Set the resistances based on the head loss type.
        if ref[:options]["headloss"] == "h-w"
            ref[:alpha] = 1.852
            ref[:resistance] = calc_resistances_hw(ref[:links])
            ref[:resistance_cost] = calc_resistance_costs_hw(ref[:links])
        elseif ref[:options]["headloss"] == "d-w"
            ref[:alpha] = 2.0
            viscosity = ref[:options]["viscosity"]
            ref[:resistance] = calc_resistances_dw(ref[:links], viscosity)
            ref[:resistance_cost] = calc_resistance_costs_dw(ref[:links], viscosity)
        else
            Memento.error(LOGGER, "Head loss formulation type is not recognized.")
        end
    end

    return refs
end
