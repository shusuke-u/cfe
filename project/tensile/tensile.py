''' field '''
n_dim = 3

''' output '''
step_out = 200

CFL_factor = 0.01

''' material '''
index = 0
eng = 1
dens = 1

''' spring '''
set_alpha = True
# spring_to_sph_ratio : 1e-2 -> 1e2, n=10
spring_to_sph_ratio = 1

kernel_support_ratio = 3

''' construction '''
# mode = "make"
# n_ptcl = 1e5
# dens = 1
# domain = [[0,0,0],[1,1,1]]
# structure = "Cartesian"
# structure = "triangle"

mode = 'read'
path = '../../../relax/result/0414/dens 3d n=1e4 h=3 thin/final.tsv'

''' tensile '''
direction = 'x'
# msa_converge_bound = 1e-4
min_step_bound = 10000
tensile_length_ratio = 0.005
poisson_converge_bound = 1e-6
side_width_ratio = 3 # ratio to mean particle spacing

dissipate= 1
# clear_vel_interval = 0
# clear_vel_lower = 1000