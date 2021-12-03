function correct_des_pipes!(data::Dict{String,<:Any})
    # Get the WaterModels portion of the data dictionary.
    wm_data = get_wm_data(data)

    # Get nondimensionalization entries, if they exist.
    base_length = get(wm_data, "base_length", 1.0)
    base_mass = get(wm_data, "base_mass", 1.0)
    base_time = get(wm_data, "base_time", 1.0)

    # Get the global head loss and viscosity parameters.
    head_loss = wm_data["head_loss"]
    viscosity = wm_data["viscosity"]

    # Define a function using the above parameters that operates on `x`.
    func! =
        x -> _correct_des_pipes!(x, head_loss, viscosity, base_length, base_mass, base_time)

    # Apply the above function to all subnetworks.
    apply_wm!(func!, data; apply_to_subnetworks = true)
end


function _correct_des_pipes!(
    data::Dict{String,<:Any},
    head_loss::String,
    viscosity::Float64,
    base_length::Float64,
    base_mass::Float64,
    base_time::Float64,
)
    capacity = calc_capacity_max(data)

    for des_pipe in values(data["des_pipe"])
        # Get common connecting node data for later use.
        node_fr = data["node"][string(des_pipe["node_fr"])]
        node_to = data["node"][string(des_pipe["node_to"])]

        # Correct various pipe properties. The sequence is important, here.
        _correct_flow_direction!(des_pipe)
        _correct_des_pipe_flow_bounds!(
            des_pipe,
            node_fr,
            node_to,
            head_loss,
            viscosity,
            capacity,
            base_length,
            base_mass,
            base_time,
        )
    end
end


function _correct_des_pipe_flow_bounds!(
    pipe::Dict{String,<:Any},
    node_fr::Dict{String,<:Any},
    node_to::Dict{String,<:Any},
    form::String,
    viscosity::Float64,
    capacity::Float64,
    base_length::Float64,
    base_mass::Float64,
    base_time::Float64,
)
    # Calculate flow bounds, which depend on the head loss formulation.
    flow_min = _calc_pipe_flow_min(
        pipe,
        node_fr,
        node_to,
        form,
        viscosity,
        capacity,
        base_length,
        base_mass,
        base_time,
    )

    flow_max = _calc_pipe_flow_max(
        pipe,
        node_fr,
        node_to,
        form,
        viscosity,
        capacity,
        base_length,
        base_mass,
        base_time,
    )

    pipe["flow_min"] = min(0.0, flow_min)
    pipe["flow_max"] = max(0.0, flow_max)

    pipe["flow_min_forward"] = max(flow_min, get(pipe, "flow_min_forward", 0.0))
    pipe["flow_max_reverse"] = min(flow_max, get(pipe, "flow_max_reverse", 0.0))

    if get(pipe, "y_min", 0.0) == 1.0
        pipe["flow_min"] = max(0.0, pipe["flow_min"])
        pipe["flow_min_forward"] = max(0.0, pipe["flow_min_forward"])
        pipe["flow_max_reverse"] = 0.0
        pipe["flow_max"] = max(0.0, pipe["flow_max"])
    elseif get(pipe, "y_max", 1.0) == 0.0
        pipe["flow_min"] = min(0.0, pipe["flow_min"])
        pipe["flow_min_forward"] = 0.0
        pipe["flow_max_reverse"] = min(0.0, pipe["flow_max_reverse"])
        pipe["flow_max"] = min(0.0, pipe["flow_max"])
    end

    pipe["flow_max"] = max(pipe["flow_min"], pipe["flow_max"])
    pipe["flow_min"] = min(pipe["flow_min"], pipe["flow_max"])
    pipe["flow_min_forward"] = max(pipe["flow_min_forward"], pipe["flow_min"])
    pipe["flow_min_forward"] = min(pipe["flow_min_forward"], pipe["flow_max"])
    pipe["flow_max_reverse"] = min(pipe["flow_max_reverse"], pipe["flow_max"])
    pipe["flow_max_reverse"] = max(pipe["flow_max_reverse"], pipe["flow_min"])

    if get(pipe, "z_max", 1.0) == 0.0
        pipe["flow_min"] = 0.0
        pipe["flow_max"] = 0.0
        pipe["flow_min_forward"] = 0.0
        pipe["flow_max_reverse"] = 0.0
    end

    @assert pipe["flow_min"] <= pipe["flow_max"]
    @assert get(pipe, "flow_min_forward", 0.0) <= max(0.0, pipe["flow_max"])
    @assert min(0.0, pipe["flow_min"]) <= get(pipe, "flow_max_reverse", 0.0)
end
