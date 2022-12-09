#input data file
junctions_mat = [
1 'o'
2 'o'
3 'o'
4 'o'
5 'o'
6 'o'
7 'r'
8 'r'
9 'r'
10 'r'
11 'r'
12 'r'
];

#edge_id type from_node to_node diameter length friction factor thermal_loss_coefficient(gamma)
pipes_mat = [
1 'o' 1 2 0.3048 1000 0.01 10
2 'o' 2 3 0.3048 1000 0.01 10
3 'o' 2 4 0.3048 1000 0.01 10
4 'o' 3 5 0.3048 1000 0.01 10
5 'o' 4 5 0.3048 1000 0.01 10
6 'o' 5 6 0.3048 1000 0.01 10
7 'r' 7 8 0.1016 1000 0.001 1e-5
8 'r' 8 9 0.1016 1000 0.001 1e-5
9 'r' 8 10 0.1016 1000 0.001 1e-5
10 'r' 9 11 0.1016 1000 0.001 1e-5
11 'r' 10 11 0.1016 1000 0.001 1e-5
12 'r' 11 12 0.1016 1000 0.001 1e-5
];

# edge_id fr to tset pset
plants_mat = [
13 12 1 #t #p
];

#edge_id fr to #q
loads_mat = [
14 2 11 3.425e6
15 3 10 3.425e6
16 4 9 3.425e6
17 5 8 3.425e6
18 6 7 3.425e6
];
