# SEM2D_READ_FAULT reads fault outputs from SEM2DPACK
#
# SYNTAX	sem2d_read_fault(model_name)
#
# INPUT		model_name	The name of directory that contains the results from SEM2DPACK
#
#               fault_name      [Flt*] 	prefix of header and data files (fault_name_sem2d.*) 
#				The default is the first FltXX_sem2d.* found 
#				in the current directory.
#
# OUTPUT	data.nx		number of fault nodes
#		data.nt		number of time samples
#		data.dt		time step
#		data.x,data.z	coordinates of fault nodes
#		data.d		slip [nx,nt]
#		data.v		slip rate
#		data.st		shear stress
#		data.sn		normal stress
#		data.mu		friction coefficient
#		data.st0	initial value of shear stress
#		data.sn0	initial value of normal stress
#		data.mu0	initial value of friction coefficient
# 
#		If output on each side of the fault (osides=T):
#  		data.d1t	displacement on side 1, fault parallel component
#  		data.d2t	displacement on side 2, fault parallel component
#  		data.d1n	displacement on side 1, fault normal component
#  		data.d2n	displacement on side 2, fault normal component
#  		data.v1t	velocity on side 1, fault parallel component
#  		data.v2t	velocity on side 2, fault parallel component
#  		data.v1n	velocity on side 1, fault normal component
#  		data.v2n	velocity on side 2, fault normal component
#
# NOTE		Fault normal components on each side of the fault are exported only in P-SV

import os
import numpy as np

def sem2d_read_fault(model_name,fault_name):
    
    # length of the tag at the begining and end of a binary record
    # in number of single precision words (4*bytes)
    LENTAG = 2; # gfortran older versions
    LENTAG = 1;
    
    # assumes header file name is FltXX_sem2d.hdr
    if not os.path.isdir(model_name):
        print("Wrong path to the model directory...")
        exit()
    headfile_exist = os.path.isfile(model_name+"/"+fault_name+"_sem2d.hdr")
    initfile_exist = os.path.isfile(model_name+"/"+fault_name+"_init_sem2d.tab")
    datafile_exist = os.path.isfile(model_name+"/"+fault_name+"_sem2d.dat")
    if (not headfile_exist):
        print("Miss head file in this directory...")
        exit()
    elif (not initfile_exist):
        print("Miss init file in this directory...")
        exit()
    elif (not datafile_exist):
        print("Miss fault data files in this directory...")
        exit()
    
    data = {}
    
    f = open(model_name+"/"+fault_name+"_sem2d.hdr")
    lines = f.readlines()
    data['nx'] = int(lines[1].split()[0])
    ndat       = int(lines[1].split()[1])
    data['nt'] = int(lines[1].split()[2])
    data['dt'] = float(lines[1].split()[3])
    xyz = []
    for line in lines[4::]:
        xyz.append(line.split())
    xyz = np.asarray(xyz).astype(np.float)
    data['x'] = xyz[:,0]
    data['z'] = xyz[:,1]

    # Read initial fault data
    f = open(model_name+"/"+fault_name+"_init_sem2d.tab")
    lines = f.readlines()
    xyz = []
    for line in lines[4::]:
        xyz.append(line.split())
    xyz = np.asarray(xyz).astype(np.float)
    data['st0'] = xyz[:,0]
    data['sn0'] = xyz[:,1]
    data['mu0'] = xyz[:,2]
    
    # Read fault data in a big matrix
    f   = open(model_name+"/"+fault_name+"_sem2d.dat", "rb")
    dt  = np.dtype((np.float32, data['nx']+2*LENTAG))
    raw = np.fromfile(f, dtype=dt)

    raw = np.reshape(raw[:,LENTAG:LENTAG+data['nx']],(int(raw.shape[0]/ndat),ndat, data['nx']));

    # Reformat each field [nx,nt]
    data['d']  = raw[:,0,:] 
    data['v']  = raw[:,1,:] 
    data['st'] = raw[:,2,:] 
    data['sn'] = raw[:,3,:] 
    data['mu'] = raw[:,4,:] 
    if (ndat == 5+4):
        data['d1t'] = raw[:,5,:] 
        data['d2t'] = raw[:,6,:] 
        data['v1t'] = raw[:,7,:] 
        data['v2t'] = raw[:,8,:] 
    elif (ndat == 5+4*2):
        data['d1t'] = raw[:,5,:] 
        data['d1n'] = raw[:,6,:] 
        data['d2t'] = raw[:,7,:] 
        data['d2n'] = raw[:,8,:] 
        data['v1t'] = raw[:,9,:] 
        data['v1n'] = raw[:,10,:] 
        data['v2t'] = raw[:,11,:] 
        data['v2n'] = raw[:,12,:] 

    return data


