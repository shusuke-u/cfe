''' field '''
n_dim = 3

''' output '''
print = True
step_out = 500

CFL_factor = 1e-3

''' material '''
index = 0
eng = 1

''' spring '''
spring_dens_ratio : 1e0 -> 1e2, n=10
# spring_dens_ratio = 1
spring_const = [5, 1, spring_dens_ratio]

spring_length_ratio = 3
kernel_support_ratio = 3

''' construction '''
# mode = "make"
# n_ptcl = 1e5
# dens = 1
# domain = [[0,0,0],[1,1,1]]
# structure = "Cartesian"
# structure = "triangle"

mode = 'read'
path = '../../../relax/result/dens 3d n=1e4 sprlen=2/final.tsv'
# path = '../../../relax/result/struct 3d n=1e4/final.tsv'
# path = '../../../../resource/distribution/3d/n=1e4/box/0.tsv'

''' tensile '''
direction = 'x'
# msa_converge_bound = 1e-4
tensile_length_ratio = 0.005
poisson_converge_bound = 1e-6
side_width_ratio = 1.5 # ratio to mean particle spacing

dissipation = 0
clear_vel_interval = 5
clear_vel_lower = 0