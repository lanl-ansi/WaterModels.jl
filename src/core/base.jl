export GenericWaterModel, build_generic_model, ids, nws

""
abstract type AbstractWaterFormulation end

"""
```
type GenericWaterModel{T<:AbstractWaterFormulation}
    model::JuMP.Model
    data::Dict{String,Any}
    setting::Dict{String,Any}
    solution::Dict{String,Any}
    var::Dict{Symbol,Any} # model variable lookup
    ref::Dict{Symbol,Any} # reference data
    ext::Dict{Symbol,Any} # user extentions
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
    model::Model

    data::Dict{String, Any}
    setting::Dict{String, Any}
    solution::Dict{String, Any}

    ref::Dict{Symbol, Any} # data reference data
    var::Dict{Symbol, Any} # JuMP variables
    con::Dict{Symbol, Any} # JuMP constraint references
    cnw::Int # current network index value

    # Extensions should define a type to hold information particular to
    # their functionality, and store an instance of the type in this
    # dictionary keyed on an extension-specific symbol.
    ext::Dict{Symbol, Any}
end

function GenericWaterModel(data::Dict{String, Any}, T::DataType;
                           ext = Dict{String, Any}(),
                           setting = Dict{String,Any}(),
                           solver = JuMP.UnsetSolver())
    ref = build_ref(data)
    var = Dict{Symbol, Any}(:nw => Dict{Int, Any}())
    con = Dict{Symbol, Any}(:nw => Dict{Int, Any}())
    cnw = minimum(keys(ref[:nw]))

    for nw_id in keys(ref[:nw])
        var[:nw][nw_id] = Dict{Symbol, Any}()
        con[:nw][nw_id] = Dict{Symbol, Any}()
    end

    wm = GenericWaterModel{T}(Model(solver = solver), data, setting,
                              Dict{String, Any}(), ref, var, con, cnw, ext)

    return wm
end

function build_generic_model(data_::Dict{String, Any}, model_constructor, post_method; kwargs...)
    data = deepcopy(data_)
    wm = model_constructor(data; kwargs...)
    post_method(wm)
    return wm
end

function build_generic_model(path::String, model_constructor, post_method; kwargs...)
    data = WaterModels.parse_file(path)
    return build_generic_model(data, model_constructor, post_method; kwargs...)
end

function build_generic_model(path::String, modification_path::String, model_constructor, post_method; kwargs...)
    data = WaterModels.parse_file(path)
    modifications = WaterModels.parse_file(modification_path)
    InfrastructureModels.update_data!(data, modifications)
    return build_generic_model(data, model_constructor, post_method; kwargs...)
end

function build_generic_model(data::Dict{String, Any}, modification_path::String, model_constructor, post_method; kwargs...)
    modifications = WaterModels.parse_file(modification_path)
    InfrastructureModels.update_data!(data, modifications)
    return build_generic_model(data, model_constructor, post_method; kwargs...)
end

function run_generic_model(data::Dict, model_constructor, solver, post_method; solution_builder = get_solution, kwargs...)
    wm = build_generic_model(data, model_constructor, post_method; kwargs...)
    return solve_generic_model(wm, solver; solution_builder = solution_builder)
end

function run_generic_model(path::String, model_constructor, solver, post_method; solution_builder = get_solution, kwargs...)
    wm = build_generic_model(path, model_constructor, post_method; kwargs...)
    return solve_generic_model(wm, solver; solution_builder = solution_builder)
end

function run_generic_model(path::String, modification_path::String, model_constructor, solver, post_method; solution_builder = get_solution, kwargs...)
    wm = build_generic_model(path, modification_path, model_constructor, post_method; kwargs...)
    return solve_generic_model(wm, solver; solution_builder = solution_builder)
end

function solve_generic_model(wm::GenericWaterModel, solver; solution_builder = get_solution)
    JuMP.setsolver(wm.model, solver)
    status, solve_time = solve(wm)
    return build_solution(wm, status, solve_time; solution_builder = solution_builder)
end

function JuMP.setsolver(wm::GenericWaterModel, solver::MathProgBase.AbstractMathProgSolver)
    JuMP.setsolver(wm.model, solver)
end

function JuMP.solve(wm::GenericWaterModel)
    status, solve_time, solve_bytes_alloc, sec_in_gc = @timed solve(wm.model, relaxation = false, suppress_warnings = true)

    try
        solve_time = getsolvetime(wm.model)
    catch
        #warn(LOGGER, "There was an issue with getsolvetime() on the solver, falling back on @timed.  This is not a rigorous timing value.");
    end

    return status, solve_time
end

"""
Returns a dict that stores commonly-used, precomputed data from of the data
dictionary, primarily for converting data types, filtering out deactivated
components, and storing system-wide values that need to be computed globally.

Some of the common keys include:

* `:pipes` -- the set of pipes in the network,
* `:reservoirs` -- the set of reservoirs in the network,
* `:junctions` -- the set of junctions in the network,
* `:valves` -- the set of valves in the network,
* `:tanks` -- the set of tanks in the network
"""
function build_ref(data::Dict{String, Any})
    refs = Dict{Symbol, Any}()
    nws = refs[:nw] = Dict{Int, Any}()

    if InfrastructureModels.ismultinetwork(data)
        nws_data = data["nw"]
    else
        nws_data = Dict{String, Any}("0" => data)
    end

    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        ref = nws[nw_id] = Dict{Symbol, Any}()

        for (key, item) in nw_data
            if isa(item, Dict{String, Any})
                item_lookup = Dict{Int, Any}([(parse(Int, k), v) for (k, v) in item])
                ref[Symbol(key)] = item_lookup
            else
                ref[Symbol(key)] = item
            end
        end

        ref[:connection] = ref[:pipes]
        ref[:ne_pipe] = filter(is_ne_pipe, ref[:pipes])
        ref[:connection_known_direction] = filter(has_known_flow_direction, ref[:connection])
        ref[:connection_unknown_direction] = filter(!has_known_flow_direction, ref[:connection])
        ref[:resistance] = calc_resistances_hw(ref[:connection])

        junction_ids = [collect(keys(ref[:junctions])); collect(keys(ref[:reservoirs]))]
        ref[:junction_connections] = Dict(i => [] for i in junction_ids)

        for (idx, connection) in ref[:connection]
            i = parse(Int, connection["node1"])
            j = parse(Int, connection["node2"])
            push!(ref[:junction_connections][i], idx)
            push!(ref[:junction_connections][j], idx)
        end
    end

    return refs
end

ids(wm::GenericWaterModel, key::Symbol) = ids(wm, wm.cnw, key)
ids(wm::GenericWaterModel, n::Int, key::Symbol) = keys(wm.ref[:nw][n][key])
ismultinetwork(wm::GenericWaterModel) = length(wm.ref[:nw]) > 1
nws(wm::GenericWaterModel) = keys(wm.ref[:nw])

# Aliases in preparation for migration to future versions of JuMP.
set_start_value = setvalue
