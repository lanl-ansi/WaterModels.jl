# # Global
R_s =  461.5 #J/kg.K
c_p_s = 1996 #specific_heat_capacity for steam
c_v_w = 4186 #specific_heat_capacity for water
c_L = 2230*1000 #latent_heat_capacity
rho_s = 0.5 #density of steam
rho_w = 1000 #density of water

#input data file
junctions_mat = [
1 'o'
2 'o'
3 'r'
4 'r'
5 'o'
6 'r'
7 'o'
8 'r'
];

#edge_id type from_node to_node diameter length friction factor thermal_loss_coefficient(gamma)
pipes_mat = [
2 'o' 1 2 0.3048 1000 0.01 10
4 'r' 3 4 0.1016 1000 0.01 10
5 'o' 2 5 0.3048 1000 0.01 10
7 'r' 6 3 0.1016 1000 0.01 10
8 'o' 5 7 0.3048 1000 0.01 10
10 'r' 8 6 0.1016 1000 0.01 10
];

# edge_id fr to tset pset
plants_mat = [
1 4 1 #t #p
]

#edge_id fr to q
loads_mat = [
3 2 3 #q
6 5 6
9 7 8
]
