# pyright: ignore

""" field """
n_dim = 3

# step_end = 0
# output_interval = step_end / 100
output_interval = 1000

periodic = [True,True,True]

CFL_factor = 0.1
kernel_support_ratio = 3

""" material """
index = 0
eng = 1
dens = 1
vel = [1,0,0] * 0.0
dissipation = 0

""" arrangement """
mode = "make"
structure = "Cartesian"
# structure = "densest"
domain = [[0,0,0],[1,1,1]]
n_ptcl = 10**(5)

""" spring """
# set_alpha = True
# spring_to_sph_ratio : 1e-2 -> 1e2, n=10
set_both = True
bulk_sound_sq = 1
# spring_const = 1 * (n_ptcl * 1e-4)**(2/n_dim-2)
spring_const = 0

""" perturbation """
amplitude = [1,0,0] * 1e-3
wave_num = 2 * pi * [1,0,0] * 1 # type: ignore
perturb_pos = True
perturb_vel = False
n_period = 5

side_width_ratio = 3