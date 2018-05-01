function parse_file(file::String; import_all = false)
    if endswith(file, ".inp")
        network_data = WaterModels.parse_epanet_file(file)
    end

    return network_data
end
