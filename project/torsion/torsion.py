''' field '''
n_dim = 3

''' output '''
step_out = 500

CFL_factor = 1e-3

''' material '''
index = 0
eng = 1

''' spring '''
sph_to_spr_ratio : 1e-1 -> 1e1, n=10
spring_const = 1

spring_length_ratio = 3
kernel_support_ratio = 3

''' construction '''
mode = "make"
n_ptcl = 1e4
dens = 1
# structure = "Cartesian"
structure = "densest"

''' torsion '''
angle = pi/12 #type: ignore
radius = 0.1
height = 1

# msa_converge_bound = 1e-3
torque_converge_bound = 1e-2

dissipation = 0
clear_vel_interval = 5
clear_vel_lower = 0