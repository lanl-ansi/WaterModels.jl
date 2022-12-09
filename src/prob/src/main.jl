# main
using JuMP, Ipopt, KNITRO
nlp_optimizer = Ipopt.Optimizer
# nlp_optimizer = KNITRO.Optimizer
model = Model(nlp_optimizer)

data_from_csv = 0;
include("parameters/global_carrier_parameters.jl")
include("parameters/default_pipe_parameters.jl")

if (data_from_csv == 1)
    include("reading_csv.jl")
    csv_junc_file = "../Kodiak_thermal_model/model_K_junction.csv"
    csv_edge_file = "../Kodiak_thermal_model/model_K_edge.csv"
    csv_plant_load_file = "../Kodiak_thermal_model/model_K_load.csv"
    csv_data = parse_csv(csv_junc_file, csv_edge_file, csv_plant_load_file)

    model_name = split(csv_junc_file,"/")[2]
    include("create_jl_file.jl")
    create_jl_file(csv_data, default_pipe_data, model_name)

    include(string("network_files/input_",model_name,".jl"))
else
    # include("network_files/input_single.jl")
    # include("network_files/input_single_double_load.jl")
    # include("network_files/input_double.jl")
    # include("network_files/input_triple.jl")`
    include("network_files/input_6.jl")
end

include("reading_network_data.jl")
include("parameters/changeable_parameters.jl")
include("reorganizing_read_data.jl")

println("***Starting Optimization***")
include("variables.jl")
include("constraints.jl")
include("objective.jl")
include("optimize.jl")
