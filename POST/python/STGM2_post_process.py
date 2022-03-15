# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 16:37:53 2022

Programme post-process-fault-2.py
    FUNCTION get_process_path()
        !!!toujours active!!!
        effectue via clavier la recherche du modele numerique
        retourne le path de la source
    FUNCTION new_model_path(local_path,mod)
        !!!toujours active!!!
        ...
    FUNCTION (path,local_path,mod)
        !!!toujours active!!!
        ...
    FUNCTION post_process_fault(process_path,fault_name):
        !!!toujours active!!!
        ...
    FUNCTION export_data_snapshot((process_path,nt,nx,nodes_x,v,st)
        !!!toujours active!!!
        ...
    FUNCTION export_data_asv()
        !!!toujours active!!!
        ...
@author: Florian Llorens / HuiHui Weng
"""


from sem2d_read_fault import sem2d_read_fault

import numpy as np
from scipy.interpolate import griddata
import scipy.ndimage.filters as filters
from os import listdir , remove , mkdir
from os.path import isfile, join, exists,isdir
import shutil
import glob

#------------------------------------------------------------------------------------
#                                '''''''path'''''''
#------------------------------------------------------------------------------------
local_path="/u/moana/user/llorens/STGM2/output/"#dossier out  
out_path="/u/moana/user/llorens/STGM2/scripts/"#dossier in
local_process_path="/u/moana/user/llorens/STGM2/scripts/newPar/"#dossier in


#------------------------------------------------------------------------------------
#                                '''''''saisie clavier du repertoire'''''''
#------------------------------------------------------------------------------------

def get_process_path():
    global local_process_path#="D:/STAGEM2/scripts/data/"#dossiers de la source : contient un seul modele de simulation
    global out_path#="D:/STAGEM2/output/"
    fichiers = [f for f in listdir(local_process_path) if not isfile(join(local_process_path, f))]
    m=[];ch=-1;i=0
    for fi in fichiers:
        i=i+1;m.append(fi);print('\n',i,'>',fi)
    while ch not in list(range(1,1+len(fichiers))):
        ch=input('>');
        try:ch=int(ch)
        except ValueError:print('!!')
    _path=out_path+m[ch-1]+'/'
    if exists(_path)==False:
        mkdir(_path);#print('mkdir : ',_path)
        
        
    return m[ch-1]
  
  
def new_model_path(mod):
    global local_path
    new_path=local_path+'new'
    if not isdir(new_path):mkdir(new_path);#print('mkdir :',new_path)
    if not isdir(new_path+'/'+mod):mkdir(new_path+'/'+mod);#print('mkdir :',new_path+'/'+mod)
    return new_path+'/'+mod

def copy_files(_path,mod):
    src=mod+'/'
    fault_name="Flt01";#print('copy : ',_path)
    shutil.copyfile(src+'Par.inp', _path+'/savePar.inp')
    shutil.copyfile(src+'info', _path+'/''info')
    shutil.copyfile(src+fault_name+'_sem2d.hdr', _path+'/'+fault_name+'_sem2d.hdr')
    shutil.copyfile(src+fault_name+'_sem2d.dat', _path+'/'+fault_name+'_sem2d.dat')
    shutil.copyfile(src+fault_name+'_init_sem2d.tab', _path+'/'+fault_name+'_init_sem2d.tab')
    shutil.copyfile(src+'MemoryInfo_sem2d.txt', _path+'/'+'MemoryInfo_sem2d.txt')

def del_files(_path):
    fichiers = [f for f in listdir(_path) if isfile(join(_path, f))]
    for fic in fichiers:
        if fic!='Par.inp':remove(_path+'/'+fic)#;print('suppr : ',fic)
    return 0

def post_process_fault(mod,fault_name):
    global local_path
    global local_process_path
    _path2=local_process_path+mod
    data = sem2d_read_fault(_path2,fault_name)   
    #V0 = float(model_name.rstrip().split("_")[12])/float(model_name.rstrip().split("_")[14])
    #V_critical = V0 * 10.0
    V_critical = 0.001
    Vs = 3330.0
    grid_size = 100.0
    
    nx       = data['nx']
    nt       = data['nt']
    dt       = data['dt']
    nodes_x  = data['x']
    v        = data['v']#vitesse de rupture ? de glissement ?
    d        = data['d']# ?
    st       = data['st']# ?
    sn       = data['sn']# ?
    mu       = data['mu']#coeff de friction ? De Lame ?
    st0      = data['st0']
    sn0      = data['sn0']
    mu0      = data['mu0']
    Rup_t    = np.zeros((nx))#temps de rupture (quand  v = v_critical )
    Rup_id   = np.zeros((nx))#indice de temps de rupture
    Heal_t   = np.zeros((nx))
    Heal_id  = np.zeros((nx))
    Vmax     = np.zeros((nx)) #vitesse 
    Gc       = np.zeros((nx)) #energie de rupture
    G0       = np.zeros((nx)) #taux de liberation d'energie statique
    Slip     = np.zeros((nx)) #glissement sur la falle
    Dtau     = np.zeros((nx)) #variation de la contrainte tangantielle
    
    #######################DATA PROCESSING####################################################################
    
    if(nt > v.shape[0]+1):
        print("Data steps is less than nt...\n Model is:", nt,v.shape[0])
        nt = v.shape[0]+1
    
    #Get rupture time for each grid
    for x in range(nx):
        for i in range(nt-2):
            if(v[i,x]<V_critical and v[i+1,x]>=V_critical):
                Rup_t[x]  = i * dt
                Rup_id[x] = i
                #print('vcrit!! :  ',Rup_t[x])
                break    
    for x in range(nx):
        Heal_t[x]  =  nt * dt
        Heal_id[x] =  nt - 1
        for i in range(nt-2):
            if(i>Rup_id[x] and v[i,x] < V_critical/1e10):
                Heal_t[x]  = i * dt
                Heal_id[x] = i
                #print('vcrit2!! :',Heal_t[x])
                break   
    for x in range(nx):# Get energies for each grid
        for i in range(nt-2):
            if(i < Rup_id[x] or i > Heal_id[x] or i==0):
                continue
            sliprate = v[i,x]
            tau      = st[i,x]
            Gc[x] = Gc[x] + dt * sliprate * tau
            if(sliprate > Vmax[x]):
                Vmax[x] = sliprate
        Gc[x] = Gc[x] - st[-1,x] * d[-1,x]
        G0[x] = 0.5 * d[-1,x] * (st[0,x] - st[-1,x])
        Slip[x] = d[-1,x]
        Dtau[x] = st[0,x] - st[-1,x]
        
    X_lower  = np.min(nodes_x)
    X_upper  = np.max(nodes_x)
    X_dim    = int((X_upper-X_lower)/grid_size+1)
    grid_x = np.mgrid[X_lower:X_upper:X_dim*1j]
    Init_t0  = griddata(nodes_x, Rup_t, grid_x, method='linear')
    grid_hl  = griddata(nodes_x,Heal_t, grid_x, method='linear')
    grid_Gc  = griddata(nodes_x,    Gc, grid_x, method='linear')
    grid_G0  = griddata(nodes_x,    G0, grid_x, method='linear')
    grid_Vmax= griddata(nodes_x,  Vmax, grid_x, method='linear')
    grid_Slip= griddata(nodes_x,  Slip, grid_x, method='linear')
    grid_Dtau= griddata(nodes_x,  Dtau, grid_x, method='linear')
    vr       = np.zeros((X_dim))
    acce     = np.zeros((X_dim))
    for i in range(X_dim-3):
        if(Init_t0[i]==0.0 or i<2): continue
        if (Init_t0[i+2]+Init_t0[i+1]-Init_t0[i-1]-Init_t0[i-2] != 0.0):
            #vr[i] =  (6*grid_size) / (Init_t0[i+2]+Init_t0[i+1]-Init_t0[i-1]-Init_t0[i-2])
            vr[i] =  (2*grid_size) / (Init_t0[i+1]-Init_t0[i-1])    
    vr = filters.gaussian_filter1d(vr,  6, mode='wrap')
    grid_Vmax = filters.gaussian_filter1d(grid_Vmax,  12, mode='wrap')   
    for i in range(X_dim):
        if(i<10 or i>X_dim-12): continue
        acce[i] =  (vr[i+1] - vr[i-1]) / (2*grid_size) * vr[i]        
    _PATH=new_model_path(mod)
    copy_files(_PATH,_path2)
    export_data_asv(_PATH,X_dim,grid_x,Init_t0,grid_hl,vr,acce,grid_Gc,grid_G0,grid_Vmax,grid_Slip,grid_Dtau)
    export_data_snapshot(_PATH,nt,nx,nodes_x,v,st)  
     
def export_data_asv(process_path,X_dim,grid_x,Init_t0,grid_hl,vr,acce,grid_Gc,grid_G0,grid_Vmax,grid_Slip,grid_Dtau):
    output= open(process_path+"/along_strike_values.dat","w")
    for i in range(X_dim):
        output.writelines(str(grid_x[i]/1e3))
        output.writelines("  ")
        output.writelines(str(Init_t0[i]))
        output.writelines("  ")
        output.writelines(str(grid_hl[i]))
        output.writelines("  ")
        output.writelines(str(vr[i]))
        output.writelines("  ")
        output.writelines(str(acce[i]))
        output.writelines("  ")
        output.writelines(str(grid_Gc[i]))
        output.writelines("  ")
        output.writelines(str(grid_G0[i]))
        output.writelines("  ")
        output.writelines(str(grid_Vmax[i]))
        output.writelines("  ")
        output.writelines(str(grid_Slip[i]))
        output.writelines("  ")
        output.writelines(str(grid_Dtau[i]))
        output.writelines("\n")
    output.close()

def export_data_snapshot(process_path,nt,nx,nodes_x,v,st):
    snap1 = int(1*nt/4.0)
    snap2 = int(2*nt/4.0)
    snap3 = int(3*nt/4.0)
    snap4=  int(4*nt/4.0)   
    output= open(process_path+"/snapshot_1.dat","w")
    for i in range(nx):
        output.writelines(str(nodes_x[i]/1e3))
        output.writelines("  ")
        output.writelines(str(v[snap1,i]))
        output.writelines("  ")
        output.writelines(str(st[snap1,i] - st[0,i]))
        output.writelines("\n")
    output.close()    
    output= open(process_path+"/snapshot_2.dat","w")
    for i in range(nx):
        output.writelines(str(nodes_x[i]/1e3))
        output.writelines("  ")
        output.writelines(str(v[snap2,i]))
        output.writelines("  ")
        output.writelines(str(st[snap2,i] - st[0,i]))
        output.writelines("\n")
    output.close()    
    output= open(process_path+"/snapshot_3.dat","w")
    for i in range(nx):
        output.writelines(str(nodes_x[i]/1e3)) 
        output.writelines("  ")
        output.writelines(str(v[snap3,i]))
        output.writelines("  ")
        output.writelines(str(st[snap3,i] - st[0,i]))
        output.writelines("\n")
    output.close()    

#------------------------------------------------------------------------------------
#                                '''programme '''
#------------------------------------------------------------------------------------    
#
#
#
# mod=get_process_path()#choix du process dans local_source_path
# #
# post_process_fault(mod,fault_name="Flt01")
# #!!!!!
# #!
# del_files(local_process_path+mod)#efface des fichier dans le local_source_path









