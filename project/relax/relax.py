''' field '''
n_dim = 3

''' output '''
step_out = 10
step_end = 500
periodic = [True, True, True]
periodic_domain = [[0,0,0],[1,1,1]]

CFL_factor = 0.1

''' material '''
index = 0
dens = 1
eng = 1
spring_to_sph_ratio = 1

kernel_support_ratio = 3

# dissipation = 1e-1
clear_vel_interval = 100
clear_vel_lower = 0
max_msa_bound = 1e-10

repulsion_factor = 1

''' construction '''
mode = 'make'
n_ptcl = 1e4
domain = [[0,0,0],[1,1,1]]
structure = 'Cartesian'
# structure = 'densest'
random_shift = 1
# structure = 'stratified'