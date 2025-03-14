#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Purpose : python interface for launching sem2dpack numerical experiments
# Hugo Lestrelin - hugo.lestrelin@hotmail.com

import numpy as np
import sys,os
import pickle,gzip
import os.path
import math
from os import path
import subprocess as sp

dir='' #working directory
cn='SF1_svw_1h15-s-5--35-aw-5' #codename of the exp
cn_1=['anglewave']#['f0','std'] #subcodename of the main parameter alterated

# verify if dir of output exist, else, create them
if not path.exists(dir+cn):
    print('Creating '+cn+' experiment.')
    sp.run(['mkdir',dir+cn])
else:
    print('Experiment '+cn+' already exists.')
if not path.exists(dir+cn+'/'+cn_1[0]):
    print('Creating '+cn_1[0]+' sub-experiment.')
    sp.run(['mkdir',dir+cn+'/'+cn_1[0]])
else:
    print('Sub-experiment '+cn_1[0]+' already exists.')

#Specific parameters 
RSss=True #initial theta proportional to Dc and Vstar 
frictionknown=False # if False, friction proportional to initial stress

#Experiment parameters 
ipar=10 #nbr of iteration

#Gaussian*sinus fct:
def gaussxsin(ampli,f0,std,onset,t):
    #sinarg=[math.radians(f0*np.pi*2*(x-onset)) for x in t]
    sinarg=f0*np.pi*2*(t-onset)
    return ampli*np.exp(-1/2*((t-onset)/std)**2)*np.sin(sinarg)

nrjref=[]
#ampliref=[]
lenexp=1000
tt=np.linspace(0,30,lenexp)
f0a=[1]#np.linspace(1.5,3,59)
for i in range(1): #59
    nrjref2=[]
    ampliref=np.linspace(0.01*9.81/(2*np.pi*f0a[i]),0.3*9.81/(2*np.pi*f0a[i]),ipar)
    for j in range(len(ampliref)):
        dd=gaussxsin(ampliref[j],f0a[i],3,15,tt)
        nrjref2.append(np.sum(dd**2)/lenexp)
    nrjref.append(nrjref2)
    #plt.plot(ampliref,nrjref2,'k',lw=4) #np.linspace(0.01,0.3,lenexp)
#plt.show()

# #ampliref=[]
# lenexp=100000
# tt=np.linspace(0,30,lenexp)
# f0a=[9.5]#np.linspace(1.5,3,59)
# for i in range(1): #59
#     nrjref2bis=[]
#     ampliref2=np.linspace(0.001*9.81/(2*np.pi*f0a[i]),3*9.81/(2*np.pi*f0a[i]),lenexp)
#     for j in range(len(ampliref2)):
#         dd=gaussxsin(ampliref2[j],f0a[i],3,15,tt)
#         nrjref2bis.append(np.sum(dd**2)/lenexp)

# newamplis=[]
# for i in range(ipar):
#     ibis=np.where(nrjref2[i]<=nrjref2bis)[0][0]
#     if np.abs(nrjref2bis[ibis]-nrjref2[i])<np.abs(nrjref2bis[ibis-1]-nrjref2[i]):
#         newamplis.append(ampliref2[ibis])
#     else:
#         newamplis.append(ampliref2[ibis-1])

# physical parameters
# fault geometric parameters
par={'angle':np.linspace(math.radians(-5),math.radians(-35),ipar) } #[math.radians(-5)]*ipar       #np.linspace(math.radians(2),math.radians(40),ipar)                                                         ### deg ; pendage faille 
par['cosangle']=[math.cos(par['angle'][i]) for i in range(ipar)]
par['sinangle']=[math.sin(par['angle'][i]) for i in range(ipar)]
par['h']=[15]*ipar#np.linspace(5,50,ipar)#[15]*ipar#np.linspace(1,100,ipar)#[20]*ipar
par['h'][i]+1-4*(par['h'][i]+1)/32                                                   ### m ; hauteur de la masse de sédiment

#global physical parameters
par['g']=[9.81]*ipar                                                                                   ### m/s² ; accélération gravitationnelle
par['poisson']=[0.2]*ipar #0.35
par['cs1']=[300]*ipar #np.linspace(300,1200,ipar)#[300,600,1200]*ipar                                  ### m.s-1 ; S-wave celerity over the fault
par['cp1']=[par['cs1'][i]*np.sqrt((1-par['poisson'][i])/(0.5-par['poisson'][i])) for i in range(ipar)] ### m.s-1 ; P-wave celerity over the fault
par['cs2']=[300]*ipar                                                                                  ### m.s-1 ; S-wave celerity below the fault
par['cp2']=[par['cs2'][i]*np.sqrt((1-par['poisson'][i])/(0.5-par['poisson'][i])) for i in range(ipar)] ### m.s-1 ; P-wave celerity below the fault
par['rho1']=[2000]*ipar #np.linspace(2000,3500,ipar)                                                   ### kg.m-3 ; density of the sediments over the fault
par['rho2']=par['rho1']#[2650]*ipar                                                                    ### kg.m-3 ; density of the sediments over the fault
par['etaKV1']=[0.1]*ipar                                                                               ### damping parameter, avoiding numerical instability
par['etaKV2']=par['etaKV1']

#source parameters
stf=True 
superpo=False
if stf:
    stfname=['synth6.5_1.tab']*ipar#['bord_1.89_15m'+'_'+str(j)+'.tab' for j in range(ipar)]#['stfpga'+str(j)+'.tab' for j in range(ipar)]
    #pore pressure param
    par['ampli']=[0.45]*ipar
    par['onset']=[15]*ipar   
    par['f0']=[1.89]*ipar
    par['std']=[3.75]*ipar
    par['fac']=[0]#np.linspace(518,841.75,ipar)#518]*ipar
    par['pga']=[0.18]#np.linspace(0.05*9.81,0.3*9.81,ipar)#[1]*ipar
else:
    par['f0']=[6]*ipar#f0a*ipar#np.geomspace(0.1,30,ipar)#[1]*ipar #300 #np.arange(0.01,30,30/ipar)            ### Hz ; ???
    par['nrj']=nrjref[0]
    par['pga']=np.linspace(0.01,0.5,ipar)#np.linspace(0.01,0.5,ipar)#[False]#np.linspace(0.01,0.3,ipar) # % of g [0.1]*ipar#
    if par['pga'][0]!=False:
        par['ampli']=[par['pga'][i]*par['g'][i]/(2*np.pi*par['f0'][i]) for i in range(ipar)]
    else:
        par['ampli']=[0.23463174999999997]#np.geomspace(0.025,2.5,ipar)#newamplis#[0.25]*ipar#np.geomspace(1*10**(-3),10,ipar) #[1]*ipar #np.linspace(1,30,ipar)              ### m.s-1 ; max amplitude of the wave at gaussian peak
    par['std']=[6]*ipar#[6]*ipar#np.linspace(8,100,ipar)#np.geomspace(0.1,100,ipar)                            ### standard deviation of the gaussian (lenght of the wave)
    par['onset']=[30]*ipar                                                                            ### sec ; when does the ampli max is acheived
    if superpo:
        tt=np.linspace(0,60,1000)
        for l in range(ipar):
            hcs=par['h'][l]*2/par['cs1'][l]
            stfref=gaussxsin(par['ampli'][l],par['f0'][l],par['std'][l],par['onset'][l],tt)
            superpo1 = list(np.zeros(np.int(hcs/(tt[1]-tt[0])+1)))
            superpo1.extend(list(stfref))
            superpo2 = list(stfref)
            superpo2.extend(list(np.zeros(np.int(hcs/(tt[1]-tt[0])+1))))
            superpotot=[]
            [superpotot.append(superpo2[i]-superpo1[i]) for i in range(len(superpo2))]
            par['ampli'][l]=par['ampli'][l]*(par['pga'][l]*par['g'][l]/(np.max(superpotot)*np.pi*2*par['f0'][l]))
    stfname=['xxxxxxx']*ipar
par['anglewave']=np.linspace(5,35,ipar)-5#np.linspace(-30,30,ipar) (0,90,ipar)   #[0]*ipar                                                                          ### deg ; angle of the incident wave

#R&S parameters
par['Dc']=[0.0004]*ipar   #0.0004                                                                             ### m ; nucleation lenght
par['Tn0']=[-(a*b*(c+1-4*(c+1)/32))*d for a,b,c,d in zip(par['rho1'],par['g'],par['h'],par['cosangle'])] #+30*1000*10   ### Pa ; initial normal stress (lithostatic) #+30*1000*10 [-592038.5090759244]*ipar#
par['Tt0']=[-(a*b*(c+1-4*(c+1)/32))*d for a,b,c,d in zip(par['rho1'],par['g'],par['h'],par['sinangle'])] #+30*1000*10  ### Pa ; initial tengantial stress (lithostatic) #+30*1000*10 [51796.65791493325]*ipar#
par['A']=[0.01]*ipar#[0.01]*ipar#np.geomspace(0.001,0.01,ipar)
afactor=np.linspace(0.5,2,ipar)  
par['B']=[0.015]*ipar#[par['A'][i]*afactor[i] for i in range(ipar)]#[0.015]*ipar#np.geomspace(0.0015,0.15,ipar)#np.arange(0.001,0.1,0.1/ipar)#0.01
par['V0']=[10**-9]*ipar#np.geomspace(1*10**-11,1*10**-4,ipar)                                          ### m.s-1 ; initial celerity of the fault
if RSss:
    if frictionknown:
        par['MuS']=[0.3]*ipar                                                                       ### inital friction
        par['Vstar']=[np.exp((c-a/b)*(1/d))*e for a,b,c,d,e in zip(par['Tt0'],par['Tn0'],par['MuS'],par['A'],par['V0'])] ### m.s-1 ; AVEC CONSIDÉRATION QUE V0 CONNU MAIS PAS VSTAR / WHAT ABOUT L'INVERSE ???
    else:
        par['Vstar']=par['V0']#[5*10**(-11)]*ipar                                                      ### m.s-1 ; AVEC CONSIDÉRATION QUE V0 CONNU MAIS PAS VSTAR / WHAT ABOUT L'INVERSE ???
        par['MuS']=[np.abs(a/b) for a,b in zip(par['Tt0'],par['Tn0'])]#[0.086]*ipar #                  ### inital friction
    par['theta']=[a/b for a,b in zip(par['Dc'],par['Vstar'])]                                          ### state variable
#PRÉVOIR ÉVOL NON RSss ????

#Visualisation parameters
par['receivers']=[3]*ipar
par['reicx1']=[0]*ipar
par['reicz1']=[-1]*ipar
par['reicx2']=[0]*ipar
par['reicz2']=[15]*ipar
par['champ']=['V']*ipar #DVS
par['itdss']=[100000000]*ipar

#Modelisation parameters
par['nelemx']=[4]*ipar #4 #16
par['nelemy']=[32]*ipar #32 # 128 #96
par['ezflt']=[4]*ipar
#par['zlim2']=[par['h'][i]-1+par['nelemx'][i]* for i in range(ipar)]
par['dx']=[2/4]*ipar
par['tempstot']=[250]*ipar #[2.5]*ipar #2.5 #60                                                         ### sec ; total lenght of experiment
par['dt']=[5*10**-6]*ipar                                                                               ### sec ; timestep
par['CFL']=[0.55]*ipar#[a*b/c for a,b,c in zip(par['dt'],par['cs1'],par['dx'])] #0.55 [0.55]*ipar#      ### Fredriech-Courant param

#checking if exp with same param already done
exprepo=dir+cn+'/'+cn_1[0]+'/'+str("{:.2e}".format(par[cn_1[0]][0]))+'___'+str("{:.2e}".format(par[cn_1[0]][-1]))
if not path.exists(exprepo):
    sp.run(['mkdir',exprepo])
else:
    print('Watch out, this experiment already exist, try another namecode / parameters !')
    exit

pickle.dump(par,open(exprepo+'/par.pkl','wb'))

for i in range(ipar):
    if i<10: cnbr='000'+str(i)
    elif i>=10 and i<100: cnbr='00'+str(i)
    elif i>=100 and i<1000: cnbr='0'+str(i)
    else: cnbr=str(i)
    sp.run(['mkdir',exprepo+'/'+cnbr]) # cn_1[j]
    if stf:
        sp.run(['cp',stfname[i],exprepo+'/'+cnbr])
    lines=[
        "SMTD evaluated as RS friction fault / parameter evaluated = "+cn, #title
        #----- Some general parameters ----------------
        "&GENERAL iexec=1, ngll=5, fmax=3.d0 , ndof=1 ,",
        "title = ' "+cn+" ', verbose='0000' , ItInfo =1000000 /",
        #----- Build the mesh ---------------------------
        "&MESH_DEF  method = 'CARTESIAN'/",
        "&MESH_CART xlim=-1,1, zlim=-1,"+str(par['h'][i])+", nelem="+str(par['nelemx'][i])+","+str(par['nelemy'][i])+", ezflt=4/", #+par['h'][i]/16
        "&MESH_CART_DOMAIN tag=1, ex= 1,"+str(par['nelemx'][i])+", ez=1,"+str(int(par['nelemy'][i]/16+3))+" /", #+0
        #"&MESH_CART_DOMAIN tag=2, ex= 1,"+str(par['nelemx'][i])+", ez="+str(int(par['nelemy'][i]/16+1))+","+str(int(par['nelemy'][i]/16+2))+" /",
        #"&MESH_CART_DOMAIN tag=3, ex= 1,"+str(par['nelemx'][i])+", ez="+str(int(par['nelemy'][i]/16+3))+","+str(int(par['nelemy'][i]/16+4))+" /",
        "&MESH_CART_DOMAIN tag=2, ex= 1,"+str(par['nelemx'][i])+", ez="+str(int(par['nelemy'][i]/16+4))+","+str(int(par['nelemy'][i]))+" /", #+5
        #---- Material parameters --------------
        "&MATERIAL tag=1, kind='ELAST' /",
        #"&MATERIAL tag=2, kind='ELAST', 'KV' /",
        "&MAT_ELASTIC rho="+str(par['rho1'][i])+", cp="+str(par['cp1'][i])+", cs="+str(par['cs1'][i])+" /",
        #"&MAT_KV eta="+str(par['etaKV1'][i])+" /",
        #"&MATERIAL tag=3, kind='ELAST', 'KV' /",
        "&MATERIAL tag=2, kind='ELAST' /", #4
        "&MAT_ELASTIC rho="+str(par['rho2'][i])+", cp="+str(par['cp2'][i])+", cs="+str(par['cs2'][i])+" /",
        #"&MAT_KV eta="+str(par['etaKV2'][i])+" /",
        #----- Boundary conditions ---------------------
        "&BC_DEF  tag = 1, kind = 'ABSORB' /", ## ADD 3 FOR UP IN CASE OF INF PLANE ; BUT DELETE FOR MODEL WITH FREE SURFACE
        #"&BC_DEF  tag = 3, kind = 'ABSORB' /",
        "&BC_DEF  tags = 2,4, kind = 'PERIOD' /",
        "&BC_DEF  tags = 5,6,   kind = 'DYNFLT' /",
        "&BC_DYNFLT friction='RSF', Tn="+str(par['Tn0'][i])+", Tt="+str(par['Tt0'][i])+", V="+str(par['V0'][i])+", pp='GAUSSXSIN', sk=0, fa=-5, cs=200, gg=100000000 /", # otd=0.00001, oxi(1)=2 oxi(2)=3, oxi(3)=1 #5000000000 #sk=0.3
        "&STF_GAUSSXSIN ampli="+str(par['ampli'][i])+" onset="+str(par['onset'][i])+" f0="+str(par['f0'][i])+" std="+str(par['std'][i])+" /",
        "&BC_DYNFLT_RSF kind=2, Dc="+str(par['Dc'][i])+", MuS="+str(par['MuS'][i])+", a="+str(par['A'][i])+", b="+str(par['B'][i])+", Vstar="+str(par['Vstar'][i])+", theta="+str(par['theta'][i])+" /",
        #---- Time scheme settings ----------------------
        "&TIME  TotalTime="+str(par['tempstot'][i])+", Courant="+str(par['CFL'][i])+", kind='newmark' /", #Dt="+str(par['dt'][i])+" /",
        #---- Sources -----------------------------------
        #---- Sources - GAUSSXSIN -----------------------
        #"&SRC_DEF stf='GAUSSXSIN', mechanism='WAVE', coord=0,0 /", #GAUSSXSIN
        "&STF_GAUSSXSIN ampli="+str(par['ampli'][i])+" onset="+str(par['onset'][i])+" f0="+str(par['f0'][i])+" std="+str(par['std'][i])+" /",
        #---- Sources - FICHIER -------------------------
        "&SRC_DEF stf='TAB', mechanism='WAVE', file='"+stfname[i]+"' /",
        "&STF_TAB file='"+stfname[i]+"' /", #
        "&SRC_WAVE angle="+str(par['anglewave'][i])+" phase='S' /",
        #----- Receivers ---------------------------------
        "&REC_LINE number = "+str(par['receivers'][i])+", first = "+str(par['reicx1'][i])+","+str(par['reicz1'][i])+", last = "+str(par['reicx2'][i])+","+str(par['reicz2'][i])+", isamp=1, field='"+par['champ'][i]+"' /",
        #--------- Plots settings ----------------------
        "&SNAP_DEF itd="+str(par['itdss'][i])+", fields ='V', components='xyz',bin=T,ps=F /"
        ]
    f=open(exprepo+'/'+cnbr+'/Par.inp','w')
    f.write('\n'.join(lines))
    f.close()
    sp.Popen("sem2dsolvepp", cwd=exprepo+'/'+cnbr+"/")


