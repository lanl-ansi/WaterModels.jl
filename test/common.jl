function build_mn_data(base_data; replicates::Int=2)
    mn_data = WaterModels.parse_file(base_data)
    return WaterModels.replicate(mn_data, replicates)
end
