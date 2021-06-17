# for my Mac
CUBIT_DIR = '/path_to_cubit_app'
curdir = 'path_to_current_directory'

###
import math, os, sys
import numpy as np
sys.path.append(CUBIT_DIR)
sys.path.append(GEOCUBIT_DIR)
import cubit
###



###
elem_size = 20




# change to current directory
command = 'cd \''+ curdir+ '\' '
cubit.cmd(command)


# by default it is 4. I avoid using all procs of my mac!
# cmd='processors 3'
# cubit.cmd(cmd)

# reset all
cubit.cmd('reset')


####################  geometry ####################


# create a 2D block (surface)
cmd = 'create surface rectangle width 1800 height 400 zplane'
cubit.cmd(cmd)


# translate to our coordinates
cmd = 'move Surface 1 x 0 y 200 z 0 include_merged'
cubit.cmd(cmd)



# ####################  mesh ####################

# set size
cmd = 'surface 1 size '+ str(elem_size)
cubit.cmd(cmd)

cmd = 'mesh surface 1'
cubit.cmd(cmd)


# refining the fault vicinity to 5 m 
cubit.cmd('refine curve 1 depth 1 size 5 bias 1.0')


# ####################  match materials and blaocks ####################
#
cubit.cmd( 'block 1 surface 1 ')
#
cubit.cmd( 'create material name \'MAT1 \' ')
#
cubit.cmd( 'block 1 material \'MAT1 \' ')

# ####################  Boundary conditions ####################
cubit.cmd('sideset 1 curve 3')
cubit.cmd('sideset 2 curve 4')
cubit.cmd('sideset 3 curve 1')
cubit.cmd('sideset 4 curve 2')

# ####################  write out in abaqus ####################
cmd = 'export abaqus "path_to_where_to_save_mesh/example_box_mesh.inp"  overwrite  everything '
cubit.cmd(cmd)
###################################################
