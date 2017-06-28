#########################################################################
#                                                                       #
# This file provides functions for interfacing with EPANET data files   #
#                                                                       #
#########################################################################

#************************************************************
#Input : .inp file from epanet                              #
#Output: Dictionary in the desired format for WaterModels.jl#
#************************************************************
function parse_epanet(file_string)
   data_string = readstring(open(file_string))
   epdata = parse_epanet_data(data_string)
   data = change_data_organization(epdata)
   return data
end


#*****************************************************************************************
#This function changes the data organization (eg: say pipes, pumps, valves combined etc.)#
#*****************************************************************************************
function change_data_organization(epdata)
   data = Dict{AbstractString,Any}()
   #mapping junctions with corresponding coordinates
   #******************************************************************************************************************
   #Note : string(symbol) can be used to convert a symbol to string (useful when index is a string instead of number"#
   #******************************************************************************************************************
   junctions = []
   id = 1
   
   for junc in epdata["junction"]
      junc_data = Dict{AbstractString,Any}()
      for coord in epdata["coordinate"]
         if junc["index"]==coord["node"]
            junc_data["name"] = junc["index"]
            junc_data["index"] = id
            junc_data["x_coord"] = coord["x-coord"]
            junc_data["y_coord"] = coord["y-coord"]
            junc_data["elevation"] = junc["elevation"]
            junc_data["demand"] = junc["demand"]
         end
      end
      push!(junctions,junc_data)                                 
#      println(junc["index"])
   id+=1

   end
   data["junction"] = junctions

   connections = []
   id = 1
   for pipe_rows in epdata["pipe"]
      connect_data = Dict{AbstractString,Any}(
            "index" => id,
            #"index" => pipe_rows["index"],
            "type" => "pipe",
            "f_junction" => pipe_rows["f_junction"],
            "t_junction" => pipe_rows["t_junction"],
            "length" => pipe_rows["length"],
            "diameter" => pipe_rows["diameter"],
            "roughness" => pipe_rows["roughness"],
            "minorloss" => pipe_rows["minorloss"],
            "status" => pipe_rows["status"],
           )
      push!(connections,connect_data)
      id+=1;
   end
   for pump_rows in epdata["pump"]
      connect_data = Dict{AbstractString,Any}(
            "index" => id,
            #"index" => pump_rows["index"],
            "type" => "pump",
            "f_junction" => pump_rows["f_junction"],
            "t_junction" => pump_rows["t_junction"],
            "parameters" => pump_rows["parameters"]
             )
      push!(connections,connect_data)
      id+=1;
   end

   for valve_rows in epdata["valve"]
      connect_data = Dict{AbstractString,Any}(
            "index" => id,
            #"index" => valve_rows["index"],
            "type" => "valve",
            "f_junction" => valve_rows["f_junction"],
            "t_junction" => valve_rows["t_junction"],
            "diameter" => valve_rows["diameter"],
            "type" => valve_rows["type"],
            "setting" => valve_rows["setting"],
            "minorloss" => valve_rows["minorloss"],
           )
      push!(connections,connect_data)
      id+=1;
   end
   data["connection"] = connections


   return data
end


#*****************************************************************************************
#Parses data strings from EPANET .inp files#
#*****************************************************************************************

function parse_epanet_data(data_string)

   data_lines = split(data_string, '\n')
   case = Dict{AbstractString,Any}(
      ##
      ##
   )


   parsed_blocks = []  
   last_index = length(data_lines)
   index = 1
   while index <= last_index
      line = strip(data_lines[index])

      if length(line) <= 0
         index = index + 1
         continue
      end

      if contains(line,"TITLE")
         index = index + 1
         line = strip(data_lines[index])
         case["title"] = strip(line)
         index = index+1 

      end

      if contains(line,"JUNCTIONS")
         regular = 1 #lines_contain semicolon at the end
         matrix_block,index = parse_blocks(line,data_lines,index,regular)
         push!(parsed_blocks,matrix_block)
      end

      if contains(line,"RESERVOIRS")
         regular=1
         matrix_block,index = parse_blocks(line,data_lines,index,regular)
         push!(parsed_blocks,matrix_block)
     end
      
     
      if contains(line,"TANKS")
         regular=1
         matrix_block,index = parse_blocks(line,data_lines,index,regular)
         push!(parsed_blocks,matrix_block)
     end

     
      if contains(line,"PIPES")
         regular=1
         matrix_block,index = parse_blocks(line,data_lines,index,regular)
         push!(parsed_blocks,matrix_block)
     end


     
      if contains(line,"PUMPS")
         regular=1
         matrix_block,index = parse_blocks(line,data_lines,index,regular)
         push!(parsed_blocks,matrix_block)
     end

     
      if contains(line,"VALVES")
         regular=1
         matrix_block,index = parse_blocks(line,data_lines,index,regular)
         push!(parsed_blocks,matrix_block)    
      end
      
      if contains(line,"COORDINATES")
         regular=0
         matrix_block,index = parse_blocks(line,data_lines,index,regular)
         #println(matrix_block)
         push!(parsed_blocks,matrix_block)
      end
      
      index += 1
   end
   #println(parsed_blocks )
   for parsed_block in parsed_blocks
      #
      #
      if parsed_block["name"] == "junctions"
         junctions = []

         for junc_row in parsed_block["data"]
            #println(junc_row)
            if length(junc_row)!=0
               
               junc_data = Dict{AbstractString,Any}(
                     "index" => parse(Int, junc_row[1]),
                     "elevation" => parse(Float64, junc_row[2]),
                     "demand" => parse(Float64, junc_row[3])
                    )
            end
            if length(junc_row)>3
               junc_data["pattern"] = parse(Float64,junc_row[4])
            end
            
            push!(junctions,junc_data)
         end
      
         case["junction"] = junctions
         #
      elseif parsed_block["name"] == "reservoirs"
         reservoirs = []


         for reserv_row in parsed_block["data"]
            if length(reserv_row)!=0
               reserv_data = Dict{AbstractString,Any}(
                     "index" => parse(Int, reserv_row[1]),
                     "head" => parse(Float64, reserv_row[2])
                     )
            end
            if length(reserv_row)>2
               reserv_data["pattern"] = parse(Float64,reserve_row[3])
            end

            push!(reservoirs,reserv_data)
         end

         case["reservoir"] = reservoirs

      elseif parsed_block["name"] == "tanks"
         tanks = []

         for tank_row in parsed_block["data"]
            if length(tank_row)!=0
               #
               tank_data = Dict{AbstractString,Any}()

               #
            end

            push!(tanks,tank_data)
         end

         case["tank"] = tanks
         
      
     
      elseif parsed_block["name"] == "pipes"
         pipes = []

         for pipe_row in parsed_block["data"]
            if length(pipe_row)!=0
               pipe_data = Dict{AbstractString,Any}(
                     "index" => parse(Int,pipe_row[1]),
                     "f_junction" => parse(Int,pipe_row[2]),
                     "t_junction" => parse(Int,pipe_row[3]),
                     "length" => parse(Float64,pipe_row[4]),
                     "diameter" => parse(Float64,pipe_row[5]),
                     "roughness" => parse(Float64,pipe_row[6]),
                     "minorloss" => parse(Float64,pipe_row[7]),
                     "status" => parse(pipe_row[8])
                    )
            end

            push!(pipes,pipe_data)
         end
         
         case["pipe"] = pipes

      elseif parsed_block["name"] == "pumps"
         pumps = []

         for pump_row in parsed_block["data"]
            if length(pump_row)!=0
               pump_data = Dict{AbstractString,Any}(
                     "index" => parse(Int,pump_row[1]),
                     "f_junction" => parse(Int,pump_row[2]),
                     "t_junction" => parse(Int,pump_row[3]),
                     "parameters" => parse(Float64,pump_row[4])
                    )                             

            end

            push!(pumps,pump_data)
         end

         case["pump"] = pumps
      elseif parsed_block["name"] == "valves"
         valves = []
         
         for valve_row in parsed_block["data"]
            if length(valve_row)! = 0
               valve_data = Dict{AbstractString,Any}(
                     "index" => parse(Int,valve_row[1]),
                     "f_junction" => parse(Int,valve_row[2]),
                     "t_junction" => parse(Int,valve_row[3]),
                     "diameter" => parse(Float64,valve_row[4]),
                     "type" => parse(Float64,valve_row[5]),
                     "setting" => parse(Float64,valve_row[6]),
                     "minorloss" => parse(Float64,valve_row[7])
                     )
            end

            push!(valves,valve_data)
         end

         case["valve"] = valves

      elseif parsed_block["name"] == "coordinates"
         coordinates = []

         for coord_row in parsed_block["data"]
            if length(coord_row)!=0
               coord_data = Dict{AbstractString,Any}(
                     "node" => parse(Int,coord_row[1]),
                     "x-coord" => parse(Float64,coord_row[2]),
                     "y-coord" => parse(Float64,coord_row[3])
                    )
                                                    
            end

            push!(coordinates,coord_data)
         end
         
         case["coordinate"] = coordinates
        
      end#elif_end
      
   end#for end

   return case 
end#func end

#*****************************************************************************************
#Parse a block (say Junctions, Pipes etc.)
#*****************************************************************************************
function parse_blocks(line,data_lines,index,regular)
   block_lines = []
   block_name = strip(replace(replace(line,'[',' '),']',' '))
   index = index+1
   column_names = split(strip(replace(data_lines[index],';',' ')))
   #println(column_names)
   index = index+1

 
   while(line!="")
      line = strip(data_lines[index])

      push!(block_lines,line)
      index +=1
   end
  #******************************************************
  #When lines don't end with ';'
   
   if (regular==0)&&(length(block_lines)!=0)
      l = size(block_lines)[1]
      block_lines1 = []
      for i in 1:l-1
         #println(block_lines[i,:])
         line1 = (join(hcat(block_lines[i,:],";")))
         push!(block_lines1,line1)
         #block_lines = join(hcat(block_lines[i,:],";"))
      end
      block_lines = block_lines1

   end
   #***************************************************
   #println(block_lines)
   block_data = join(block_lines,' ')
   block_data = replace(strip(block_data),'\t',' ')
   #println(block_data)
   block_data_rows = split(block_data,';')

   #println(block_data_rows)
   block_data_rows = block_data_rows[1:length(block_data_rows)-1]
   #println(block_data_rows)
   matrix_block = []
   empty_warning = 0

   for row in block_data_rows
      row_items = split_line(strip(row))
      if (empty_warning==0)&&(length(row_items)!=length(column_names))
      empty_warning = 1
      end
      push!(matrix_block,row_items)
   end
  
   if size(matrix_block,1)==0
      println("Warning : All fields in $block_name are empty")
   elseif empty_warning == 1
      println("Warning : Empty fields in $block_name,are shifted to the right most column and ignored")

   end

   block_dict = Dict("name" => lowercase(block_name), "data" => matrix_block, "column_names"=>column_names)
     
   index = index-1
   return block_dict,index
end


#*****************************************************************************************
#Split line into substring array
#*****************************************************************************************
single_quote_expr = r"\'((\\.|[^\'])*?)\'"

function split_line(ep_line::AbstractString)
    if ismatch(single_quote_expr, ep_line)
        # splits a string on white space while escaping text quoted with "'"
        # note that quotes will be stripped later, when data typing occurs

        #println(ep_line)
        tokens = []
        while length(ep_line) > 0 && ismatch(single_quote_expr, ep_line)
            #println(ep_line)
            m = match(single_quote_expr, ep_line)

            if m.offset > 1
                push!(tokens, ep_line[1:m.offset-1])
            end
            push!(tokens, replace(m.match, "\\'", "'")) # replace escaped quotes

            ep_line = ep_line[m.offset+length(m.match):end]
        end
        if length(ep_line) > 0
            push!(tokens, ep_line)
        end
        #println(tokens)

        items = []
        for token in tokens
            if contains(token, "'")
                push!(items, strip(token))
            else
                for parts in split(token)
                    push!(items, strip(parts))
                end
            end
        end
        #println(items)

        #return [strip(ep_line, '\'')]
        return items
    else
        return split(ep_line)
    end
end


