# pyright: ignore

""" field """
n_dim = 3

# step_end = 0
# output_interval = step_end / 100
output_interval = 10

periodic = [True,True,True]

# CFL_factor = 0.005
CFL_factor = 0.1
kernel_support_ratio = 3.1

""" material """
index = 0
eng = 1
dens = 1
vel = [1,0,0] * 0.0
dissipate = 0

""" elastic moduli """
poisson = 0.3

""" arrangement """
mode = "make"
structure = "Cartesian"
# structure = "densest"
domain = [[0,0,0],[1,1,1]]
n_ptcl = 50**n_dim

# mode = 'read'
# path = '../../../relax/result/dens 3d n=1e4 sprlen=2/final.tsv'

""" perturbation """
ampl = 1e-3
wave_num = 2 * pi * 2

n_period = 5

side_width_ratio = 2