# conf.py

# number of dimensions
n_dim = 2
# number of particles
n_ptcl = 10000
# total mass
mass = 9

''' spring force '''
att_coef = 1
rep_coef = 1
# spring length limit to bond
spring_length_limit = 1.5

kernel_support_radius = 5

''' material '''
index = 0
reuse_time = True

n_out = 50
time_end = 2

dis_coef = 0
internal_energy = 1

''' perturbation '''
shift = [ 1e-2, 0e-2, 0 ]

''' construction '''
# construct_mode = 'read'
# path to the file in which initial distribution is written
# path = '../../res/distr/2d/n=1e3/box/0.csv'
# path = '../../res/distr/3d/n=1e4/box/0.csv'
# path = '../../res/distr/3d/n=10/0.csv'

mode = 'make'
distribution = 'Cartesian'

box = [ [ -1, -1, -1 ], [ 2, 2, 2 ] ]
energy = 1

periodic_boundary = [False, True, False]
# box_for_periodic_boundary = [ [ -0.5, 0.5 ], [ -0.5, 0.5 ], [ 0, 0 ] ]
