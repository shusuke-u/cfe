''' field '''
n_dim = 3

''' output '''
step_out = 100
step_end = 1e5

CFL_factor = 0.001

''' material '''
index = 0
eng = 1

''' spring '''
# spring_dens_ratio : 1e-6 -> 1e-4, n=10
# spring_dens_ratio = 0
spring_const = [3, 1e0, 1e0]
dens_diff_std = 1e0
spring_length_ratio = 2
kernel_support_ratio = 3

''' construction '''
mode = "make"
n_ptcl = 1e4
dens = 1
domain = [[0, 0, 0],[1, 0.1, 0.5]]
structure = "Cartesian"
# structure = "triangle"

# mode = 'read'
# path = '../../../relax/result/dens 3d n=1e4 sprlen=3/final.tsv'

''' oscillation '''
init_vel = 1