## Sets
SetNO = Set([1,2])
SetNR = Set([3,4])
SetN = union(SetNO, SetNR)
PlantSet = Set([1])
SetPO = Set([2])
LoadSet = Set([3])
SetPR = Set([4])

SetN = union(SetNO, SetNR)
SetE = union(SetPO, SetPR, LoadSet, PlantSet) #currently no pumps and plant modeled as a setpoint
SetP = union(SetPO, SetPR)

# to node for an edge
N_to = zeros(Int,4)
N_to[1] = 1
N_to[2] = 2
N_to[3] = 3
N_to[4] = 4

#from nodes for an edge
N_fr = zeros(Int,4)
N_fr[1] = 4
N_fr[2] = 1
N_fr[3] = 2
N_fr[4] = 3

# edges entering a junction
E_in = []
push!(E_in,Set([1]))
push!(E_in,Set([2]))
push!(E_in,Set([3]))
push!(E_in,Set([4]))

P_in = []
push!(P_in, Set())
push!(P_in, Set([2]))
push!(P_in, Set())
push!(P_in, Set([4]))


L_in =[]
push!(L_in, Set())
push!(L_in, Set())
push!(L_in, Set([3]))
push!(L_in, Set())

PL_in = []
push!(PL_in, Set([1]))
push!(PL_in, Set())
push!(PL_in, Set())
push!(PL_in, Set())


E_out = []
push!(E_out,Set([2]))
push!(E_out,Set([3]))
push!(E_out,Set([4]))
push!(E_out,Set([1]))

P_out = []
push!(P_out, Set([2]))
push!(P_out, Set())
push!(P_out, Set([4]))
push!(P_out, Set())

L_out =[]
push!(L_out, Set())
push!(L_out, Set([3]))
push!(L_out, Set())
push!(L_out, Set())


PL_out = []
push!(PL_out, Set())
push!(PL_out, Set())
push!(PL_out, Set())
push!(PL_out, Set([1]))
