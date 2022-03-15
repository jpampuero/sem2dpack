# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 02:20:21 2022

@author: tynes
"""

import numpy as np
import matplotlib.pyplot as plt
from os import system,popen,chdir
import subprocess
import sys


from STGM2_ import get_path_STGM2
from STGM2_ import get_Par_data
from STGM2_ import convertir_Par_data,reconvertir_Par_data
from STGM2_ import export_Par_data
from STGM2_ import print_,print__
from STGM2_ import rsf_mu
from STGM2_ import steady_state
from STGM2_post_process import post_process_fault,del_files

def update_(D,ref,newValue):
    D[ref]=newValue
    D=rsf_mu(D)
    return D

def update_dist_(x_dim,x,D,ref):
    D['TtXn']=len(x)
    if ref=='theta':D['thetaXn']=len(x);D['TtXn']=len(x)
    D[ref] =[*x_dim,*x]
    if ref=='theta':
        D=rsf_mu(D)
    return D

def calcul_theta(ss,Dc,V):
    x=[0]*len(ss);i=0
    if type(Dc)!=list:#ss vecteur  Dc,V scalair
        for e in ss:x[i]=e*Dc/V ;i=i+1
    elif type(Dc)==type(ss):#ss,Dc vecteur  V scalair
        for e,f in zip(ss,Dc):x[i]=e*f/V;i=i+1
    return x

       
def new_sim(dic,temp_,ti):
    path_,out_dir=export_Par_data(temp_,dic,ti)
    chdir(path_[:-8])
    myPopen = subprocess.Popen('sem2dsolve > info', shell = True, stdout = subprocess.PIPE, encoding = 'ascii')
    lines=[]
    while True:
      line = myPopen.stdout.readline()
      lines.append(line)
      if line == '' and myPopen.poll() is not None:
          break
    returnStatus = myPopen.poll()
    if returnStatus != 0:raise RuntimeError('Problem')
    post_process_fault(out_dir,fault_name="Flt01")
    del_files(path_[:-8])#efface des fichier dans le local_source_path        

#---------------------------------------------------------------------------
#                                '''initialisation : choix du template '''
#------------------------------------------------------------------------------------
#
#
models_path=get_path_STGM2()
m=-1
template_=models_path[m]
#
#
##------------------------------------------------------------------------------------
#                                '''debut programme : exporter en masse '''
#------------------------------------------------------------------------------------
#-PARAMETRES in-----------------------------------------
#
#
print('================Template======================')
print('path :',template_)
D=get_Par_data(template_)
print_(D)
print('=============END-Template=====================')
D=convertir_Par_data(D)
#
##-PARAMETRES update-----------------------------------------
pc=5.5564446
#____________________________________test_DC_____________________________________________________________________
#Dc_=[0.1*D['Dc'],0.2*D['Dc'],0.3*D['Dc'],0.4*D['Dc'],0.5*D['Dc'],0.6*D['Dc'],0.7*D['Dc'],0.8*D['Dc'],0.9*D['Dc']]
#s_s=[2,2,1e-3]
#t_dim=[1500,10000]
#for e in Dc_:
#    D=get_Par_data(template_)
#    D=convertir_Par_data(D)
#    #####ajout des modification
#    update_(D,'Dc',e)
#    t=calcul_theta(s_s,D['Dc'],D['V']);print('theta(Dc,V) = ',t)
#    D=update_dist_(t_dim,t,D,ref='theta') 
#    D['Tt'][2]=D['Tt'][2]*(1+pc/100)#ajout de contrainte dans la premiere zone (pc%)
#    new_sim(D,template_,ti='Dc')
#__________________________________________________________________________________________________________________

#_____________________________________test_dBarr_________________________________________________________________
# s_s=[1,1,1e-3]
# #t_dim=[[1500,5000],[1500,6000],[1500,7000],[1500,8000],[1500,9000],[1500,10000],[1500,11000],[1500,12000],[1500,13000],[1500,14000]]
# t_dim=[]
# for e in (1000*np.arange(1.5,16,0.5)):
#     t_dim.append([1500,e])
# for e in t_dim:
#     D=get_Par_data(template_)
#     D=convertir_Par_data(D)
#     #####ajout des modification
#     t=calcul_theta(s_s,D['Dc'],D['V'])
#     D=update_dist_(e,t,D,ref='theta');print('s-s:',steady_state(D))
#     D['Tt'][2]=D['Tt'][2]*(1+pc/100)#ajout de contrainte dans la premiere zone (pc%)
#     new_sim(D,template_,ti='dBarr')
# #__________________________________________________________________________________________________________________

#_____________________________________test_DcBarr_________________________________________________________________
s_s=[1,1,1e-3]
Dc_=[]
#Dc_=[[3.7e-2,3.7e-2,3.7e-1],[3.7e-2,3.7e-2,3.7e-2],[3.7e-2,3.7e-2,3.7e-3],[3.7e-2,3.7e-2,3.7e-4],[3.7e-2,3.7e-2,3.7e-5],[3.7e-2,3.7e-2,3.7e-6]]
for e in (3.7e-5 * np.arange(1,10000,20)):
    Dc_.append([3.7e-2,3.7e-2,e])

t_dim=[1500,5000]
for dc in Dc_:
    D=get_Par_data(template_)
    D=convertir_Par_data(D)
    #####ajout des modification
    t=calcul_theta(s_s,dc,D['V'])
    D=update_dist_(t_dim,dc,D,ref='Dc')
    D=update_dist_(t_dim,t,D,ref='theta');print('s-s:',steady_state(D))
    D['Tt'][2]=D['Tt'][2]*(1+pc/100)#ajout de contrainte dans la premiere zone (pc%)
    new_sim(D,template_,ti='DcTheta')
#__________________________________________________________________________________________________________________

#
# #
#------------------------------------------------------------------------------------
#                                '''debut programme : exporter 1  '''
#------------------------------------------------------------------------------------
#  
#
# # -PARAMETRES in-----------------------------------------

# print('\nIN : Par.inp')
# print_(D)
# D=convertir_Par_data(D)
# #
# #-PARAMETRES updates------------------------------------
# #
# Dc=3.7e-2
# t=[97,1,1e-3]
# t_dim=[1500,10000]
# #
# #D['Dc']    = D['Dc']     *  1.0
# #
# update_(D,'Dc',Dc)
# #update_(D,'Vstar',1e-6)
# #update_(D,'Tn',-120e6)
# #update_(D,'V',1e-3)
# #update_(D,'a',0.00625)
# #update_(D,'b',0.0086)
# #update_(D,'MuS',0.6)
# #update_(D,'Vc',0)
# #
# #
# update_dist_(t_dim,t,D,ref='theta')
# #
# print("stady state : ",steady_state(D))
# #
# ##-PARAMETRES out---------------------------------------
# export_Par_data(template_,D,ti='Dc')




    



