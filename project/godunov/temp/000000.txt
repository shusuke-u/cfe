# conf.py

# number of dimensions
n_dim = 2
# number of particles
n_ptcl = 1000
# total mass
mass = 1
ptcl_mass = mass / n_ptcl


''' spring '''
att_coef = 0
rep_coef = 1e2
# spring length limit to bond
spring_length_limit = 1.5

kernel_support_radius = 3

# temporary
max_spr_coef = max( att_coef, rep_coef )
mean_spr_coef = ( att_coef * rep_coef )**0.5

''' Courant condition '''
time_interval =  2 * ( ptcl_mass / max_spr_coef )**0.5 * 1e-3

''' material '''
index = 0
reuse_time = True

step_out = 1000

time_end = 50

dis_coef = 0
internal_energy = 1

construct_mode = 'make_Cartesian'
box = [ [ 0, 0, 0 ], [ 1, 1, 1 ] ]
energy = 1

periodic_boundary = [ True, False, False ]
# box_for_periodic_boundary = [ [ -0.5, 0.5 ], [ -0.5, 0.5 ], [ 0, 0 ] ]
