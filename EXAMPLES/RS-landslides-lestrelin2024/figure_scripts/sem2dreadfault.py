#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import sys,os
import os.path
from os import path
import subprocess as sp
from os import listdir
from os.path import isfile, join
import pickle,gzip
import matplotlib.pyplot as plt


def sem2d_read_fault(model_name,fault_name):

    # length of the tag at the begining and end of a binary rec 
    # in number of single precision words (4*bytes)
    LENTAG = 2; # gfortran older versions
    LENTAG = 1;

    # assumes header file name is FltXX_sem2d.hdr
    if not os.path.isdir(model_name):
        print("Wrong path to the model directory...")
        quit()
    headfile_exist = os.path.isfile(model_name+"/"+fault_name+"_sem2d.hdr")
    initfile_exist = os.path.isfile(model_name+"/"+fault_name+"_init_sem2d.tab")
    datafile_exist = os.path.isfile(model_name+"/"+fault_name+"_sem2d.dat")
    if (not headfile_exist):
        print("Miss head file in this directory...")
        quit()
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
    data['theta0'] = xyz[:,3]

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
    data['theta'] = raw[:,5,:]
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


def chercher_VLSdiff1(exp,cb,exptype):
    vdmat=[]
    for i in range(len(cb)):
        name=exp+cb[i]+'_resultplots_'+exptype
        vdmat_data=pickle.load(open(name+'/vlsdiff.pkl','rb'))
        if i==0:
            ref_len=len(vdmat_data)
        print(name)
        print(ref_len-len(vdmat_data))
        if (ref_len-len(vdmat_data))==0: 
            vdmat.append(vdmat_data)
        else:
            vdmat.append(vdmat_data[:ref_len-len(vdmat_data)])
    return vdmat#np.array(vdmat,dtype=object)

def chercher_VLSdiff2(exp,exptype):
    vdmat=[]
    ld=os.listdir(exp+'/plots/')
    for i in range(len(ld)):
        #print(ld[i])
        if ld[i][6:6+len(exptype)] == exptype:
        #if ld[i][12:12+len(exptype)] == exptype:
            name=exp+'/plots/'+ld[i]
            #vdmat_data=pickle.load(open(name+'/vlsdiff.pkl','rb'))
            vdmat_data=pickle.load(open(name+'/norml2.pkl','rb'))
            nbr=ld[i][7+len(exptype):]
            #ld2=os.listdir(exp+'/'+exptype+'/'+str("{:.2e}".format(nbr))+'___'+str("{:.2e}".format(nbr))+'/')

            if i==0:
                ref_len=len(vdmat_data)
            if len(vdmat_data) < ref_len:
                print(name)
                print(ref_len-len(vdmat_data))
                for i in range(len(vdmat)):
                    del(vdmat[i][len(vdmat_data)-ref_len])
                ref_len=len(vdmat_data)
            if (ref_len-len(vdmat_data))==0: 
                vdmat.append(vdmat_data)
            else:
                vdmat.append(vdmat_data[:ref_len-len(vdmat_data)])
    return vdmat#np.array(vdmat,dtype=object)

def chercher_VLSdiff3(exp,exptype):
    vdmat={}
    vdmat['factor']=[]
    vdmat['data']=[]
    ld=os.listdir(exp+'/plots/')
    for i in range(len(ld)):
        #print('this is ld',ld[i])
        if ld[i][6:6+len(exptype)] == exptype:
            #print('this is ld passing the condition',ld[i])
        #if ld[i][12:12+len(exptype)] == exptype:
            name=exp+'/plots/'+ld[i]
            #vdmat_data=pickle.load(open(name+'/vlsdiff.pkl','rb'))
            vdmat_data=pickle.load(open(name+'/norml2.pkl','rb'))
            #print(len(vdmat_data))
            nbr=float(ld[i][7+len(exptype):])
            if i==0:
                vdmat['factor'].append(nbr)
                vdmat['data'].append(vdmat_data)
            else:
                j=0
                while vdmat['factor'][j]<nbr:
                    j+=1
                    #print(j)
                    if j>=len(vdmat['factor']):
                        break
                vdmat['factor'].insert(j,nbr)
                vdmat['data'].insert(j,vdmat_data)
    return vdmat,ld#np.array(vdmat,dtype=object)

def chercher_theta(exp,exptype):
    vdmat={}
    vdmat['factor']=[]
    vdmat['data']=[]
    ld=os.listdir(exp+'/'+exptype+'/')
    for i in range(len(ld)):
        name=exp+'/'+exptype+'/'+ld[i]
        nbr=float(ld[i][0:8])
        #vdmat_data=pickle.load(open(name+'/vlsdiff.pkl','rb'))
        ld2=os.listdir(name)
        vdmat_data=[]
        stfplus=[]
        facfac=[]
        for k in range(len(ld2)):
            try:  
                data=pickle.load(open(name+'/'+ld2[k]+'/data.pkl','rb'))
                if k==0:
                    facfac.append(float(ld2[k]))
                    vdmat_data.append(data['theta'][:,int(data['nx']/2)][-1]/data['theta'][:,int(data['nx']/2)][0])
                else:
                    j=0
                    while facfac[j]<float(ld2[k]):
                        j+=1
                        #print(j)
                        if j>=len(facfac):
                            break
                    facfac.insert(j,float(ld2[k]))
                    vdmat_data.insert(j,data['theta'][:,int(data['nx']/2)][-1]/data['theta'][:,int(data['nx']/2)][0])
            except(NotADirectoryError):
                #print('error')
                continue
        #print(len(vdmat_data))
        if i==0:
            vdmat['factor'].append(nbr)
            vdmat['data'].append(vdmat_data)
        else:
            j=0
            while vdmat['factor'][j]<nbr:
                j+=1
                #print(j)
                if j>=len(vdmat['factor']):
                    break
            vdmat['factor'].insert(j,nbr)
            vdmat['data'].insert(j,vdmat_data)
    return vdmat,ld#np.array(vdmat,dtype=object)

def chercher_vcommetheta(exp,exptype):
    vdmat={}
    vdmat['factor']=[]
    vdmat['data']=[]
    ld=os.listdir(exp+'/'+exptype+'/')
    for i in range(len(ld)):
        name=exp+'/'+exptype+'/'+ld[i]
        nbr=float(ld[i][0:8])
        #vdmat_data=pickle.load(open(name+'/vlsdiff.pkl','rb'))
        ld2=os.listdir(name)
        vdmat_data=[]
        stfplus=[]
        facfac=[]
        for k in range(len(ld2)):
            try:  
                data=pickle.load(open(name+'/'+ld2[k]+'/data.pkl','rb'))
                if k==0:
                    facfac.append(float(ld2[k]))
                    vdmat_data.append(data['v'][:,int(data['nx']/2)][-1]/data['v'][:,int(data['nx']/2)][0])
                else:
                    j=0
                    while facfac[j]<float(ld2[k]):
                        j+=1
                        #print(j)
                        if j>=len(facfac):
                            break
                    facfac.insert(j,float(ld2[k]))
                    vdmat_data.insert(j,data['v'][:,int(data['nx']/2)][-1]/data['v'][:,int(data['nx']/2)][0])
            except(NotADirectoryError):
                #print('error')
                continue
        #print(len(vdmat_data))
        if i==0:
            vdmat['factor'].append(nbr)
            vdmat['data'].append(vdmat_data)
        else:
            j=0
            while vdmat['factor'][j]<nbr:
                j+=1
                #print(j)
                if j>=len(vdmat['factor']):
                    break
            vdmat['factor'].insert(j,nbr)
            vdmat['data'].insert(j,vdmat_data)
    return vdmat,ld#np.array(vdmat,dtype=object)

def chercher_vmaxcommetheta(exp,exptype):
    vdmat={}
    vdmat['factor']=[]
    vdmat['data']=[]
    ld=os.listdir(exp+'/'+exptype+'/')
    for i in range(len(ld)):
        name=exp+'/'+exptype+'/'+ld[i]
        nbr=float(ld[i][0:8])
        #vdmat_data=pickle.load(open(name+'/vlsdiff.pkl','rb'))
        ld2=os.listdir(name)
        vdmat_data=[]
        stfplus=[]
        facfac=[]
        for k in range(len(ld2)):
            try:  
                data=pickle.load(open(name+'/'+ld2[k]+'/data.pkl','rb'))
                if k==0:
                    facfac.append(float(ld2[k]))
                    vdmat_data.append(np.max(data['v'][:,int(data['nx']/2)])/data['v'][:,int(data['nx']/2)][0])
                else:
                    j=0
                    while facfac[j]<float(ld2[k]):
                        j+=1
                        #print(j)
                        if j>=len(facfac):
                            break
                    facfac.insert(j,float(ld2[k]))
                    vdmat_data.insert(j,np.max(data['v'][:,int(data['nx']/2)])/data['v'][:,int(data['nx']/2)][0])
            except(NotADirectoryError):
                #print('error')
                continue
        #print(len(vdmat_data))
        if i==0:
            vdmat['factor'].append(nbr)
            vdmat['data'].append(vdmat_data)
        else:
            j=0
            while vdmat['factor'][j]<nbr:
                j+=1
                #print(j)
                if j>=len(vdmat['factor']):
                    break
            vdmat['factor'].insert(j,nbr)
            vdmat['data'].insert(j,vdmat_data)
    return vdmat,ld#np.array(vdmat,dtype=object)

def chercher_tauexp(exp,exptype):
    vdmat={}
    vdmat['factor']=[]
    vdmat['data']=[]
    ld=os.listdir(exp+'/'+exptype+'/')
    for i in range(len(ld)):
        name=exp+'/'+exptype+'/'+ld[i]
        nbr=float(ld[i][0:8])
        #vdmat_data=pickle.load(open(name+'/vlsdiff.pkl','rb'))
        ld2=os.listdir(name)
        vdmat_data=[]
        stfplus=[]
        facfac=[]
        for k in range(len(ld2)):
            try:  
                data=pickle.load(open(name+'/'+ld2[k]+'/data.pkl','rb'))
                if k==0:
                    facfac.append(float(ld2[k]))
                    vdmat_data.append(np.max(data['st'][:,int(data['nx']/2)]))
                else:
                    j=0
                    while facfac[j]<float(ld2[k]):
                        j+=1
                        #print(j)
                        if j>=len(facfac):
                            break
                    facfac.insert(j,float(ld2[k]))
                    vdmat_data.insert(j,np.max(data['st'][:,int(data['nx']/2)]))
            except(NotADirectoryError):
                #print('error')
                continue
        #print(len(vdmat_data))
        if i==0:
            vdmat['factor'].append(nbr)
            vdmat['data'].append(vdmat_data)
        else:
            j=0
            while vdmat['factor'][j]<nbr:
                j+=1
                #print(j)
                if j>=len(vdmat['factor']):
                    break
            vdmat['factor'].insert(j,nbr)
            vdmat['data'].insert(j,vdmat_data)
    return vdmat,ld#np.array(vdmat,dtype=object)

def chercher_tauinc_super(exp,exptype):
    vdmat={}
    vdmat['factor']=[]
    vdmat['data']=[]
    vdmat['xaxis']=[]
    ld=os.listdir(exp+'/'+exptype+'/')
    for i in range(len(ld)):
        #print('this is ld',ld[i])
        #if ld[i][6:6+len(exptype)] == exptype:
        #print('this is ld passing the condition',ld[i])
        #if ld[i][12:12+len(exptype)] == exptype:
        name=exp+'/'+exptype+'/'+ld[i]
        nbr=float(ld[i][0:8])
        #vdmat_data=pickle.load(open(name+'/vlsdiff.pkl','rb'))
        ld2=os.listdir(name)
        vdmat_data=[]
        vdmat_xaxis=[]
        stfplus=[]
        facfac=[]
        par=pickle.load(open(name+'/par.pkl','rb'))
        print(i,nbr)
        for k in range(len(ld2)):
            try:  
                stftot=pickle.load(open(name+'/'+ld2[k]+'/stf.pkl','rb'))
                hcs=par['h'][int(ld2[k])]*2/par['cs1'][int(ld2[k])]
                superpo1 = list(np.zeros(np.int(hcs/(stftot['t'][1]-stftot['t'][0])+1)))
                superpo1.extend(list(stftot['stf']))
                superpo2 = list(stftot['stf'])
                superpo2.extend(list(np.zeros(np.int(hcs/(stftot['t'][1]-stftot['t'][0])+1))))
                superpotot=[]
                [superpotot.append(superpo2[i]-superpo1[i]) for i in range(len(superpo2))]
                pga_inc=(2*np.pi*par['f0'][int(ld2[k])])*np.max(superpotot)
                tau_inc=par['rho1'][int(ld2[k])]*par['cs1'][int(ld2[k])]*np.max(superpotot)
                if k==0:
                    facfac.append(float(ld2[k]))
                    vdmat_data.append(tau_inc)
                    vdmat_xaxis.append(pga_inc)
                else:
                    j=0
                    while facfac[j]<float(ld2[k]):
                        j+=1
                        #print(j)
                        if j>=len(facfac):
                            break
                    facfac.insert(j,float(ld2[k]))
                    vdmat_data.insert(j,tau_inc)
                    vdmat_xaxis.insert(j,pga_inc)
            except(NotADirectoryError):
                print('?')
                continue
        #print(len(vdmat_data))
        if i==0:
            vdmat['factor'].append(nbr)
            vdmat['data'].append(vdmat_data)
            vdmat['xaxis'].append(vdmat_xaxis)
        else:
            j=0
            while vdmat['factor'][j]<nbr:
                j+=1
                #print(j)
                if j>=len(vdmat['factor']):
                    break
            vdmat['factor'].insert(j,nbr)
            vdmat['data'].insert(j,vdmat_data)
            vdmat['xaxis'].insert(j,vdmat_xaxis)
    return vdmat,ld#np.array(vdmat,dtype=object)

def chercher_tauinc(exp,exptype):
    vdmat={}
    vdmat['factor']=[]
    vdmat['data']=[]
    ld=os.listdir(exp+'/'+exptype+'/')
    for i in range(len(ld)):
        #print('this is ld',ld[i])
        #if ld[i][6:6+len(exptype)] == exptype:
        #print('this is ld passing the condition',ld[i])
        #if ld[i][12:12+len(exptype)] == exptype:
        name=exp+'/'+exptype+'/'+ld[i]
        nbr=float(ld[i][0:8])
        #vdmat_data=pickle.load(open(name+'/vlsdiff.pkl','rb'))
        ld2=os.listdir(name)
        vdmat_data=[]
        stfplus=[]
        facfac=[]
        par=pickle.load(open(name+'/par.pkl','rb'))
        for k in range(len(ld2)):
            try:  
                stftot=pickle.load(open(name+'/'+ld2[k]+'/stf.pkl','rb'))
                tau_inc=par['rho1'][int(ld2[k])]*par['cs1'][int(ld2[k])]*np.max(stftot['stf'])
                if k==0:
                    facfac.append(float(ld2[k]))
                    vdmat_data.append(tau_inc)
                else:
                    j=0
                    while facfac[j]<float(ld2[k]):
                        j+=1
                        #print(j)
                        if j>=len(facfac):
                            break
                    facfac.insert(j,float(ld2[k]))
                    vdmat_data.insert(j,tau_inc)
            except(NotADirectoryError):
                continue
        #print(len(vdmat_data))
        if i==0:
            vdmat['factor'].append(nbr)
            vdmat['data'].append(vdmat_data)
        else:
            j=0
            while vdmat['factor'][j]<nbr:
                j+=1
                #print(j)
                if j>=len(vdmat['factor']):
                    break
            vdmat['factor'].insert(j,nbr)
            vdmat['data'].insert(j,vdmat_data)
    return vdmat,ld#np.array(vdmat,dtype=object)

def chercher_pgainc_super(exp,exptype):
    vdmat={}
    vdmat['factor']=[]
    vdmat['data']=[]
    ld=os.listdir(exp+'/'+exptype+'/')
    for i in range(len(ld)):
        #print('this is ld',ld[i])
        #if ld[i][6:6+len(exptype)] == exptype:
        #print('this is ld passing the condition',ld[i])
        #if ld[i][12:12+len(exptype)] == exptype:
        name=exp+'/'+exptype+'/'+ld[i]
        nbr=float(ld[i][0:8])
        #vdmat_data=pickle.load(open(name+'/vlsdiff.pkl','rb'))
        ld2=os.listdir(name)
        vdmat_data=[]
        stfplus=[]
        facfac=[]
        par=pickle.load(open(name+'/par.pkl','rb'))
        for k in range(len(ld2)):
            try:  
                stftot=pickle.load(open(name+'/'+ld2[k]+'/stf.pkl','rb'))
                hcs=par['h'][int(ld2[k])]*2/par['cs1'][int(ld2[k])]
                superpo1 = list(np.zeros(np.int(hcs/(stftot['t'][1]-stftot['t'][0])+1)))
                superpo1.extend(list(stftot['stf']))
                superpo2 = list(stftot['stf'])
                superpo2.extend(list(np.zeros(np.int(hcs/(stftot['t'][1]-stftot['t'][0])+1))))
                superpotot=[]
                [superpotot.append(superpo2[i]-superpo1[i]) for i in range(len(superpo2))]
                pga_inc=(2*np.pi*par['f0'][int(ld2[k])])*np.max(superpotot)
                if k==0:
                    facfac.append(float(ld2[k]))
                    vdmat_data.append(pga_inc)
                else:
                    j=0
                    while facfac[j]<float(ld2[k]):
                        j+=1
                        #print(j)
                        if j>=len(facfac):
                            break
                    facfac.insert(j,float(ld2[k]))
                    vdmat_data.insert(j,pga_inc)
            except(NotADirectoryError):
                continue
        #print(len(vdmat_data))
        if i==0:
            vdmat['factor'].append(nbr)
            vdmat['data'].append(vdmat_data)
        else:
            j=0
            while vdmat['factor'][j]<nbr:
                j+=1
                #print(j)
                if j>=len(vdmat['factor']):
                    break
            vdmat['factor'].insert(j,nbr)
            vdmat['data'].insert(j,vdmat_data)
    return vdmat,ld#np.array(vdmat,dtype=object)

def chercher_stfint(exp,exptype):
    vdmat={}
    vdmat['factor']=[]
    vdmat['data']=[]
    ld=os.listdir(exp+'/'+exptype+'/')
    for i in range(len(ld)):
        #print('this is ld',ld[i])
        #if ld[i][6:6+len(exptype)] == exptype:
        #print('this is ld passing the condition',ld[i])
        #if ld[i][12:12+len(exptype)] == exptype:
        name=exp+'/'+exptype+'/'+ld[i]
        nbr=float(ld[i][0:7])
        #vdmat_data=pickle.load(open(name+'/vlsdiff.pkl','rb'))
        ld2=os.listdir(name)
        vdmat_data=[]
        stfplus=[]
        facfac=[]
        for k in range(len(ld2)):
            try:  
                stftot=pickle.load(open(name+'/'+ld2[k]+'/stf.pkl','rb'))
                stftot['stf'][np.where(stftot['stf']<0)]=0
                stf=stftot['stf']
                stfint=[stf[0]]
                [stfint.append(stfint[i]+stf[i+1]*(stftot['t'][-1]-stftot['t'][-2])) for i in range(len(stf)-1)]
                first45=np.where(stfint>0.45*np.max(stfint))
                first55=np.where(stfint>0.55*np.max(stfint))
                gradstf=(0.55*np.max(stfint)-0.45*np.max(stfint))/(stftot['t'][first55[0][0]]-stftot['t'][first45[0][0]])
                if k==0:
                    facfac.append(float(ld2[k]))
                    vdmat_data.append(gradstf)
                    stfplus.append(stfint)
                else:
                    j=0
                    while facfac[j]<float(ld2[k]):
                        j+=1
                        #print(j)
                        if j>=len(facfac):
                            break
                    facfac.insert(j,float(ld2[k]))
                    vdmat_data.insert(j,gradstf)
            except(NotADirectoryError):
                continue
        #print(len(vdmat_data))
        if i==0:
            vdmat['factor'].append(nbr)
            vdmat['data'].append(vdmat_data)
        else:
            j=0
            while vdmat['factor'][j]<nbr:
                j+=1
                #print(j)
                if j>=len(vdmat['factor']):
                    break
            vdmat['factor'].insert(j,nbr)
            vdmat['data'].insert(j,vdmat_data)
    return vdmat,ld#np.array(vdmat,dtype=object)


def chercher_param(exp,exptype,par,parindex):
    ld=os.listdir(exp+'/'+exptype)
    params=[]
    for i in range(len(ld)):
        pari=pickle.load(open(exp+'/'+exptype+'/'+ld[i]+'/par.pkl','rb'))
        params.append(pari[par][parindex])
    return params

def print_exp(exp_nbr,ou,cache):
    j=exp_number
    if j<10: cnbr='000'+str(j)
    elif j>=10 and j<100: cnbr='00'+str(j)
    elif j>=100 and j<1000: cnbr='0'+str(j)
    else: cnbr=str(j)
    file_src=ou+'/'+cnbr
    try:
        listd3=os.listdir(file_src)
        print('going in : '+file_src)
        if any(e in 'data.pkl' for e in listd3): # check output already written
            data=pickle.load(open(file_src+'/data.pkl','rb'))
        else:
            print('Creating output in : '+file_src)
            for k in range(len(listd3)):
                if listd3[k][0:3]=='Flt':
                    Fltname=listd3[k][0:5]
            data=sem2d_read_fault(file_src,Fltname)
            pickle.dump(data,open(file_src+'/data.pkl','wb'),protocol=4)
        if len(data['v'])>2:
            Vnorm=data['v']/par['V0'][j]
            stress_eff=(data['st'][:,int(data['nx']/2)]/par['Tt0'][j]+data['sn'][:,int(data['nx']/2)]/par['Tn0'][j])/2
            xt=np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache)
            
            ## fig vitesse norm
            plt.figure()
            plt.plot(xt,Vnorm[cache:,int(data['nx']/2)])
            plt.xlabel("time (s)")
            plt.ylabel("Velocity normalized")
            plt.yscale('log')
            plt.title("Velocity middle fault -- "+cn_1[0]+"="+str("{:.2e}".format(par[cn_1[0]][j])))
            plt.savefig(dir+cn+'_'+nameoutput+'/velocity/vel_'+cnbr+'.png',bbox_inches='tight',dpi=300,format='png')
            plt.close()
            
            ## fig stress tangential
            plt.figure()
            plt.plot(xt,data['st'][cache:,int(data['nx']/2)]/par['Tt0'][j])
            plt.xlabel("time (s)")
            plt.ylabel("Tangential stress normalized")
            plt.title("Tangential stress middle fault -- "+cn_1[0]+"="+str("{:.2e}".format(par[cn_1[0]][j])))
            plt.savefig(dir+cn+'_'+nameoutput+'/tengential_stress/st_'+cnbr+'.png',bbox_inches='tight',dpi=300,format='png')
            plt.close()

            ## fig comparasion solution analytique
            plt.figure()
            plt.plot(np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache),(-np.log(Vnorm[cache:,int(data['nx']/2)])*par['A'][j]),lw=2,label='-log(V/V0)*A')
            plt.plot(np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache),((data['st'][cache:,int(data['nx']/2)])/par['Tn0'][j]),linestyle='--',lw=1,label='tau/sig0')
            plt.legend()
            plt.xlabel("time (s)")
            plt.ylabel("tau/sig0 ou -Alog(V/V0)")
            plt.title("tau/sig0 vs. -Alog(V/V0) -- "+cn_1[0]+"="+str("{:.2e}".format(par[cn_1[0]][j])))
            plt.savefig(dir+cn+'_'+nameoutput+'/vs_linear-sol/VLS_'+cnbr+'.png',bbox_inches='tight',dpi=300,format='png')
            plt.close()

            ## fig soustraction solution analytique
            plt.figure()
            plt.plot(np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache),(-np.log(Vnorm[cache:,int(data['nx']/2)])*par['A'][j]-(data['st'][cache:,int(data['nx']/2)])/par['Tn0'][j]),lw=2)
            plt.legend()
            plt.xlabel("time (s)")
            plt.ylabel("-Alog(V/V0) - Tstress normalized ")
            plt.title("VLS diff-- "+cn_1[0]+"="+str("{:.2e}".format(par[cn_1[0]][j])))
            plt.savefig(dir+cn+'_'+nameoutput+'/vs_linear-sol/VLSdiff_'+cnbr+'.png',bbox_inches='tight',dpi=300,format='png')
            plt.close()

            ## 
            if len(exp_nbr)>1:
                Vnew=Vnorm[cache:,int(data['nx']/2)]
                mask=np.where(Vnew<0)[0]
                for m in range(len(mask)):
                    Vnew[mask[m]]=10**-9
                VLSmean=np.mean(np.mean((-np.log(Vnew)*par['A'][j])/((data['st'][cache:,int(data['nx']/2)])/par['Tn0'][j])))
                VLSratio.append(VLSmTsressean)
                whencheck.append(j)
                VLSdiffmean=np.mean(np.mean(-np.log(Vnorm[cache:,int(data['nx']/2)])*par['A'][j]-(data['st'][cache:,int(data['nx']/2)])/par['Tn0'][j]))
                VLSdiff.append(VLSdiffmean)
    except(ValueError):
        print('ValueError')
    return VLSdiff

def sem2d_read_STF(model_name):
    
    # length of the tag at the begining and end of a binary record
    # in number of single precision words (4*bytes)
    LENTAG = 2; # gfortran older versions
    LENTAG = 1;
    
    # assumes header file name is FltXX_sem2d.hdr
    if not os.path.isdir(model_name):
        print("Wrong path to the model directory...")
        quit()
    initfile_exist = os.path.isfile(model_name+"/SourcesTime_sem2d.tab")
    if (not initfile_exist):
        print("Miss init file in this directory...")
        exit()
    
    data = {}
    # Read initial fault data
    f = open(model_name+"/SourcesTime_sem2d.tab")
    #f = open("SourcesTime_sem2d.tab")
    lines = f.readlines()
    xyz = []
    for line in lines[4::]:
        xyz.append(line.split())
    xyz = np.asarray(xyz).astype(np.float)
    data['t'] = xyz[:,0]
    data['stf'] = xyz[:,1]
    return data