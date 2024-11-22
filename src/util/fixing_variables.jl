function fix_variables_to_relaxed_solutions(wm::AbstractWaterModel, sol_val_dict::Dict{String, Any}; nw::Int=nw_id_default)




    fixing_vector = sol_val_dict["fixing_vector"]
    fixing_all_pipes_direction = fixing_vector[1]
    fixing_pumps_direction = fixing_vector[2]
    fixing_build = fixing_vector[3]
    fixing_activation = fixing_vector[4]
    
    understanding_partitions = true

    if(fixing_all_pipes_direction == true)
        if(haskey(ref(wm,nw),:pipe))
            if(nw == 1)
                println("fixing pipe flow direction")
            end
            for (a, pipe) in ref(wm, nw, :pipe)
                val = sol_val_dict["nw"][string(nw)]["pipe"][string(a)]["y"]
                # println("fixing pipe flow direction nw $nw pipe $a to $val")
                y_pipe = var(wm, nw, :y_pipe, a)
                JuMP.unset_binary(y_pipe)
                JuMP.fix(y_pipe, val; force = true)


                # if(understanding_partitions == true)
                #     lambda_pipe = var(wm, nw, :lambda_p_pipe)
                #     x_pipe = var(wm, nw, :x_p_pipe)
                # end
            end
            if(haskey(ref(wm,nw),:short_pipe))
                if(nw == 1)
                    println("fixing short pipe flow direction")
                end
                for (a, short_pipe) in ref(wm, nw, :short_pipe)
                    val = sol_val_dict["nw"][string(nw)]["short_pipe"][string(a)]["y"]
                    val = Int(round(val))
                    y_short_pipe = var(wm, nw, :y_short_pipe, a)
                    JuMP.unset_binary(y_short_pipe)
                    JuMP.fix(y_short_pipe, val; force = true)
                end
            end
    
            if(haskey(ref(wm,nw),:ne_short_pipe))
                if(nw == 1)
                    println("fixing ne short pipe direction")
                end
                for (a, ne_short_pipe) in ref(wm, nw, :ne_short_pipe)
                    val = sol_val_dict["nw"][string(nw)]["ne_short_pipe"][string(a)]["y"]
                    val = Int(round(val))
                    y_ne_short_pipe = var(wm, nw, :y_ne_short_pipe, a)
                    JuMP.unset_binary(y_ne_short_pipe)
                    JuMP.fix(y_ne_short_pipe, val; force = true)
                end
            end
    
        end
    end


    if(fixing_pumps_direction == true)
        
        if(haskey(ref(wm,nw),:pump))
            if(nw == 1)
                println("fixing pump flow direction")
            end
            for (a, pump) in ref(wm, nw, :pump)
                val = sol_val_dict["nw"][string(nw)]["pump"][string(a)]["y"]
                val = Int(round(val))
                y_pump = var(wm, nw, :y_pump, a)
                JuMP.unset_binary(y_pump)
                JuMP.fix(y_pump, val; force = true)
            end
        end
        if(haskey(ref(wm,nw),:ne_pump))
            if(nw == 1)
                println("fixing ne pump flow direction")
            end
            for (a, ne_pump) in ref(wm, nw, :ne_pump)
                val = sol_val_dict["nw"][string(nw)]["ne_pump"][string(a)]["y"]
                val = Int(round(val))
                y_ne_pump = var(wm, nw, :y_ne_pump, a)
                JuMP.unset_binary(y_ne_pump)
                JuMP.fix(y_ne_pump, val; force = true)
            end
        end

        # if(haskey(ref(wm,nw),:valve) && haskey(sol_val_dict["nw"][string(nw)], "valve"))
        #     for (a, valve) in ref(wm, nw, :valve)
        #             yval = sol_val_dict["nw"][string(nw)]["valve"][string(a)]["y"]
        #             yval = Int(round(yval))
        #             y_valve = var(wm, nw, :y_valve, a)
        #             JuMP.unset_binary(y_valve)
        #             JuMP.fix(y_valve, yval; force = true)
        #         end
        #     end
    end

    if(fixing_build == true)
        if(haskey(ref(wm,nw),:ne_pump))
            if(nw == 1)
                println("fixing ne pump build status")
            end
            for (a, ne_pump) in ref(wm, nw, :ne_pump)
                val = sol_val_dict["nw"][string(nw)]["ne_pump"][string(a)]["build_status"]
                val = Int(round(val))
                x_ne_pump = var(wm, nw, :x_ne_pump, a)
                JuMP.unset_binary(x_ne_pump)
                JuMP.fix(x_ne_pump, val; force = true)
            end
        end
        if(haskey(ref(wm,nw),:ne_short_pipe))
            if(nw == 1)
                println("fixing ne short pipe build status")
            end
            for (a, ne_short_pipe) in ref(wm, nw, :ne_short_pipe)
                val = sol_val_dict["nw"][string(nw)]["ne_short_pipe"][string(a)]["status"]
                val = Int(round(val))
                z_ne_short_pipe = var(wm, nw, :z_ne_short_pipe, a)
                JuMP.unset_binary(z_ne_short_pipe)
                JuMP.fix(z_ne_short_pipe, val; force = true)
            end
        end
    end

    if(fixing_activation == true)
        if(haskey(ref(wm,nw),:pump))
            if(nw == 1)
                println("fixing pump activation status")
            end
            for (a, pump) in ref(wm, nw, :pump)
                val = sol_val_dict["nw"][string(nw)]["pump"][string(a)]["status"]
                val = Int(round(val))
                z_pump = var(wm, nw, :z_pump, a)
                JuMP.unset_binary(z_pump)
                JuMP.fix(z_pump, val; force = true)
            end
        end

        if(haskey(ref(wm,nw),:ne_pump))
            if(nw == 1)
                println("fixing ne pump activation status")
            end
            for (a, ne_pump) in ref(wm, nw, :ne_pump)
                val = sol_val_dict["nw"][string(nw)]["ne_pump"][string(a)]["status"]
                val = Int(round(val))
                z_ne_pump = var(wm, nw, :z_ne_pump, a)
                JuMP.unset_binary(z_ne_pump)
                JuMP.fix(z_ne_pump, val; force = true)
            end
        end
    end


    # if(haskey(ref(wm,nw),:valve) && haskey(sol_val_dict["nw"][string(nw)], "valve"))
    #     for (a, valve) in ref(wm, nw, :valve)
    #         val = sol_val_dict["nw"][string(nw)]["valve"][string(a)]["status"]
    #         val = Int(round(val))
    #         z_valve = var(wm, nw, :z_valve, a)
    #         JuMP.unset_binary(z_valve)
    #         JuMP.fix(z_valve, val; force = true)
    #     end

    # end

    
end
