function fix_variables_to_relaxed_solutions(wm::AbstractWaterModel, sol_val_dict::Dict{String, Any}; nw::Int=nw_id_default, continuous_fixing::Bool = false)
# println("Using relaxed solutions")
# println(sol_val_dict)

###############################
# Integer Variables ###########
###############################
fixing_direction = true
    if(haskey(ref(wm,nw),:pipe))
        # println("fixing pipe flow direction")
        for (a, pipe) in ref(wm, nw, :pipe)
            if(fixing_direction == true)
                val = sol_val_dict["nw"][string(nw)]["pipe"][string(a)]["y"]
                y_pipe = var(wm, nw, :y_pipe, a)
                JuMP.unset_binary(y_pipe)
                JuMP.fix(y_pipe, val; force = true)
            end
            # partition = pipe["flow_partition"]
            # qp = var(wm, nw, :qp_pipe, a)
            # qn = var(wm, nw, :qn_pipe, a)
            # qval = sol_val_dict["nw"][string(nw)]["pipe"][string(a)]["q"]
            # for i in 1:length(partition)-1
            #     if(qval >= partition[i] && qval <= partition[i+1])
            #         if(qval >=0)
            #             JuMP.set_lower_bound(qp,partition[i])
            #             JuMP.set_upper_bound(qp,partition[i+1])
            #         else
            #             JuMP.set_lower_bound(qn,partition[i])
            #             JuMP.set_upper_bound(qn,partition[i+1])
            #         end
            #     end
            # end
        end
    end

    if(haskey(ref(wm,nw),:short_pipe))
        # println("fixing pipe flow direction")
        for (a, short_pipe) in ref(wm, nw, :short_pipe)
            if(fixing_direction == true)
                val = sol_val_dict["nw"][string(nw)]["short_pipe"][string(a)]["y"]
                val = Int(round(val))
                y_short_pipe = var(wm, nw, :y_short_pipe, a)
                # println("Unsetting binary for y pipe")
                JuMP.unset_binary(y_short_pipe)
                JuMP.fix(y_short_pipe, val; force = true)
            end

        end
    end
    if(haskey(ref(wm,nw),:pump))
        println("fixing pump activation status")
        # for (a, pump) in ref(wm, nw, :pump)
        #     val = sol_val_dict["nw"][string(nw)]["pump"][string(a)]["status"]
        #     val = Int(round(val))
        #     z_pump = var(wm, nw, :z_pump, a)
        #     JuMP.unset_binary(z_pump)
        #     JuMP.fix(z_pump, val; force = true)
        # end
        println("fixing pump flow direction")
        for (a, pump) in ref(wm, nw, :pump)
            # if(fixing_direction == true)
            #     val = sol_val_dict["nw"][string(nw)]["pump"][string(a)]["y"]
            #     val = Int(round(val))
            #     y_pump = var(wm, nw, :y_pump, a)
            #     JuMP.unset_binary(y_pump)
            #     JuMP.fix(y_pump, val; force = true)
            # end

            # partition = pump["flow_partition"]
            # qp_pump = var(wm, nw, :qp_pump, a)
            # # qn_pump = var(wm, nw, :qn_pump, a)
            # qpumpval = val = sol_val_dict["nw"][string(nw)]["pump"][string(a)]["q"]
            # for i in 1:length(partition)-1
            #     if(qpumpval >= partition[i] && qpumpval <= partition[i+1])
            #         if(qpumpval >=0)
            #             JuMP.set_lower_bound(qp_pump,partition[i])
            #             JuMP.set_upper_bound(qp_pump,partition[i+1])
            #         else
            #             JuMP.set_lower_bound(qn_pump,partition[i])
            #             JuMP.set_upper_bound(qn_pump,partition[i+1])
            #         end
            #     end
            # end
        end
    end
    #
    #
    if(haskey(ref(wm,nw),:ne_pump))
        # println("fixing ne pump build status")
        for (a, ne_pump) in ref(wm, nw, :ne_pump)
            val = sol_val_dict["nw"][string(nw)]["ne_pump"][string(a)]["build_status"]
            val = Int(round(val))
            x_ne_pump = var(wm, nw, :x_ne_pump, a)
            JuMP.unset_binary(x_ne_pump)
            JuMP.fix(x_ne_pump, val; force = true)
        end
        # println("fixing ne pump activation status")
        for (a, ne_pump) in ref(wm, nw, :ne_pump)
            val = sol_val_dict["nw"][string(nw)]["ne_pump"][string(a)]["status"]
            val = Int(round(val))
            z_ne_pump = var(wm, nw, :z_ne_pump, a)
            JuMP.unset_binary(z_ne_pump)
            JuMP.fix(z_ne_pump, val; force = true)

            # # println("fixing ne pump flow direction")
            if(fixing_direction == true)
                val = sol_val_dict["nw"][string(nw)]["ne_pump"][string(a)]["y"]
                val = Int(round(val))
                y_ne_pump = var(wm, nw, :y_ne_pump, a)
                JuMP.unset_binary(y_ne_pump)
                JuMP.fix(y_ne_pump, val; force = true)
            end
        end
    end
    #
    if(haskey(ref(wm,nw),:ne_short_pipe))
        # println("fixing ne short pipe flow direction")
        for (a, ne_short_pipe) in ref(wm, nw, :ne_short_pipe)
            # if(fixing_direction == true)
            #     val = sol_val_dict["nw"][string(nw)]["ne_short_pipe"][string(a)]["y"]
            #     val = Int(round(val))
            #     y_ne_short_pipe = var(wm, nw, :y_ne_short_pipe, a)
            #     JuMP.unset_binary(y_ne_short_pipe)
            #     JuMP.fix(y_ne_short_pipe, val; force = true)
            # end

        # println("fixing ne short pipe build status")
            val = sol_val_dict["nw"][string(nw)]["ne_short_pipe"][string(a)]["status"]
            val = Int(round(val))
            z_ne_short_pipe = var(wm, nw, :z_ne_short_pipe, a)
            JuMP.unset_binary(z_ne_short_pipe)
            JuMP.fix(z_ne_short_pipe, val; force = true)

        end
    end
    if(haskey(ref(wm,nw),:valve))
        for (a, valve) in ref(wm, nw, :valve)
            val = sol_val_dict["nw"][string(nw)]["valve"][string(a)]["status"]
            val = Int(round(val))
            z_valve = var(wm, nw, :z_valve, a)
            JuMP.unset_binary(z_valve)
            JuMP.fix(z_valve, val; force = true)
        end
        for (a, valve) in ref(wm, nw, :valve)
            if(fixing_direction == true)
                yval = sol_val_dict["nw"][string(nw)]["valve"][string(a)]["y"]
                yval = Int(round(yval))
                y_valve = var(wm, nw, :y_valve, a)
                JuMP.unset_binary(y_valve)
                JuMP.fix(y_valve, yval; force = true)
            end
        end
    end
end
