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
from sem2dreadfault import sem2d_read_fault, chercher_VLSdiff1, chercher_VLSdiff2, chercher_VLSdiff3, print_exp, sem2d_read_STF, chercher_param, chercher_stfint, chercher_theta, chercher_tauexp, chercher_tauinc,chercher_tauinc_super,chercher_pgainc_super
import matplotlib as mpl
from matplotlib import colors


def pganalytichugo(f,param):
    A=[]
    for i in range(len(f)):
        x= f[i] * ( np.pi * 4 )
        z= complex(0,x)
        A.append(  (param['tau_c'] * param['cs'] / param['mu']) * np.abs(z/2 / (np.exp(-z*param['h']/ (param['cs']))-1))) #/10
        #A=   param['tau_c'] / (2 * param['h'] * param['rho'] ) * 1 / np.abs(np.exp(z)-1)
    return A


def pganalytihugo_full(f,param):
    A=[]
    for i in range(len(f)):
        x= f[i] * ( np.pi * 2 )
        z= complex(0,x)
        w0=np.log(np.pi*(param['Dc']/(param['v0']*param['duration']))**2)-np.log(np.log(np.pi*(param['Dc']/(param['v0']*param['duration']))**2)) #Lambert function
        w0=10.84 # cf. I0-1 graph
        #A.append((param['a']*param['sig0']*(1+np.abs(z)*param['theta0'])-param['b']*param['sig0']+param['mu']*param['v0']/(2*param['cs'])*(1-np.exp(-np.abs(z)*2*param['h']/param['cs']))*(1+np.abs(z)*param['theta0']))*param['cs']/param['mu'] * np.abs(z) / (1 - np.exp(-np.abs(z)*2*param['h']/param['cs'])))
        #A.append(-w0*(param['a']*param['sig0']-param['b']*param['sig0']/(1+np.abs(z)*param['theta0'])+param['mu']*param['v0']/(2*param['cs'])*(1-np.exp(-np.abs(z)*2*param['h']/param['cs'])))*param['cs']/param['mu'] * np.abs(z) / (1 - np.exp(-np.abs(z)*2*param['h']/param['cs'])))
                
        A.append(w0*np.abs((param['a']*param['sig0']-param['b']*param['sig0']/(1+z*param['theta0'])+param['mu']*param['v0']/(2*param['cs'])*(1-np.exp(-z*2*param['h']/param['cs'])))*param['cs']/param['mu'] * z / (1 - np.exp(-z*2*param['h']/param['cs']))))
        
        #A.append(np.abs((param['a']*param['sig0']*(1+z*param['theta0'])-param['b']*param['sig0']+param['mu']*param['v0']/(2*param['cs'])*(1-np.exp(-z*2*param['h']/param['cs']))*(1+z*param['theta0']))*param['cs']/param['mu'] * z / (1 - np.exp(-z*2*param['h']/param['cs']))))
        #A.append(np.abs(((1+z*param['theta0'])*(param['a']*param['sig0']+param['mu']*param['v0']/(2*param['cs'])*(1-np.exp(-(z)*2*param['h']/param['cs'])))-param['b']*param['sig0']) * (param['cs'] / param['mu'] * (z / (1 - np.exp(-z*2*param['h']/param['cs'])))))) #/10
    return A

def pganalytihugo_full_thetabis(f,param):
    A=[]
    for i in range(len(f)):
        x= f[i] * ( np.pi * 2 )
        z= complex(0,x)
        A.append(np.abs(-1/1e6*((1+z*param['theta0'])*(param['a']*param['sig0']+param['mu']*param['v0']/(2*param['cs'])*(+np.exp(-(z)*2*param['h']/param['cs'])-1))-param['b']*param['sig0']) * (param['cs'] / param['mu'] * (z / ( -np.exp(-z*2*param['h']/param['cs'])+1))))) #/10
    return A

def pganalytihugo_full_thetabis2(f,param):
    A=[]
    for i in range(len(f)):
        x= f[i] * ( np.pi * 2 )
        z= np.abs(complex(0,x))
        A.append(-1/10000000*((1+z*param['theta0'])*(param['a']*param['sig0']+param['mu']*param['v0']/(2*param['cs'])*(-np.exp(-(z)*2*param['h']/param['cs'])+1))-param['b']*param['sig0']) * (param['cs'] / param['mu'] * (z / ( -np.exp(-z*2*param['h']/param['cs'])+1)))) #/10
    return A


def pganalytijp_full_theta(f,param):
    A=[]
    for i in range(len(f)):
        x= f[i] * ( np.pi * 2 )
        z= complex(0,x)
        A.append(np.abs(((1+z*param['theta0'])*(param['a']*param['sig0']/param['v0']+param['mu']/(2*param['cs'])*(1-np.exp(-(z)*2*param['h']/param['cs'])))-param['b']*param['sig0']/param['v0']) * (param['cs'] / param['mu'] * (z / (1 - np.exp(-z*2*param['h']/param['cs'])))))) #/10
    return A


dir='/home/hugo/Documents/phd_local/sem2dpack_pptest/' #working directory
cn='svs_pga50ref' #codename of the exp
cn_1=['f0']#'B'#,'f0','std'] #subcodename of the parameter alterated
linear_sol=True
cache=0
exptype=1

### Matrix fig de l'écart à la solution analytique
#vdmat=chercher_VLSdiff1(dir+cn[:-1],list(map(str,np.arange(4,15))),str(exptype))
#vdmat=chercher_VLSdiff(dir+cn[:-1],list(map(str,np.arange(1,15))),str(exptype))
#vdmat2,ldplot2=chercher_stfint(dir+cn,cn_1[0])#chercher_VLSdiff3(dir+cn,cn_1[0])
vdmat3,ldplot3=chercher_theta(dir+cn,cn_1[0])#chercher_VLSdiff3(dir+cn,cn_1[0])
vdmat2,ldplot2=chercher_tauexp(dir+cn,cn_1[0])#chercher_VLSdiff3(dir+cn,cn_1[0])
#vdmat,ldplot=chercher_tauinc(dir+cn,cn_1[0])#chercher_VLSdiff3(dir+cn,cn_1[0])
#vdmat4,ldplot4=chercher_pgainc_super(dir+cn,cn_1[0])#chercher_VLSdiff3(dir+cn,cn_1[0])
vdmat,ldplot=chercher_tauinc_super(dir+cn,cn_1[0])#chercher_VLSdiff3(dir+cn,cn_1[0])
# maxx=np.max(vdmat['data'])
# maxreal=np.max(maxx)
# for i in range(len(vdmat['data'])):
#     vdmat['data'][i] = [x / maxreal for x in vdmat['data'][i]]

yaxis=[]
xaxis=[]
ld=os.listdir(dir+cn+'/'+cn_1[0])
for i in range(len(vdmat3['factor'])):
    for j in range(len(ld)):
        if vdmat3['factor'][i]==float(ld[j][:8]):
            #print(float(ld[j][:8]))
            par=pickle.load(open(dir+cn+'/'+cn_1[0]+'/'+ld[j]+'/par.pkl','rb'))
            yaxis.append(par['f0'])
            xaxis.append(par['pga']*10)
            # amplipga=[]
            # for k in range(len(par['f0'])):
            #     amplipga.append(par['ampli'][k]*(2*np.pi*par['f0'][k]))
            # xaxis.append(amplipga)
            #xaxis.append(par['ampli'])
    

norm = colors.TwoSlopeNorm(vmin=0, vcenter=0.5, vmax=1)#np.max(vdmat)) 0.18/0.36
norm2 = colors.TwoSlopeNorm(vmin=110000, vcenter=115424, vmax=120424)
norm4 = colors.TwoSlopeNorm(vmin=70000, vcenter=75000, vmax=80000)
norm3 = colors.LogNorm(vmin=200, vmax=500000)#np.max(vdmat)) 0.18/0.36
#norm = colors.LogNorm(vmin=np.min(vdmat['data']), vmax=np.max(vdmat['data']))#np.max(vdmat)) 0.18/0.36
#norm = colors.TwoSlopeNorm(vmin=0.01, vcenter=0.05,vmax=0.1)#np.max(vdmat)) 0.18/0.36

fig, axes = plt.subplots(nrows=1, figsize=(7, 5), constrained_layout=True)
pcm = axes.pcolormesh(xaxis, yaxis, vdmat3['data'],norm=norm, cmap=plt.cm.seismic, rasterized=True)


#pcm = axes.pcolormesh(par['pgv']*100, amplis, vdmat,norm=norm, cmap=plt.cm.seismic, rasterized=True) #vmin=-0.2, vmax=0.01 #
#pcm = axes.pcolormesh(nrjref2, amplis, vdmat,norm=norm, cmap=plt.cm.seismic, rasterized=True)
#pcm = axes.pcolormesh(vdmat,norm=norm, cmap=plt.cm.seismic, rasterized=True)
clbr=fig.colorbar(pcm, ax=axes, pad=0.01)
#clbr.ax.set_title(r"$\frac{\theta}{\theta_{0}}$",y=0.5,fontsize=16,rotation=0)
#clbr.ax.text(2.5,0.5,r"$\tau_{inc}$",fontsize=16)
clbr.ax.text(2.5,0.5,r"$\frac{\theta}{\theta_{0}}$",fontsize=16)
clbr.ax.text(0.1,0.02,"unstable",fontsize=16,rotation='vertical',color='white')
clbr.ax.text(0.1,0.84,"stable",fontsize=16,rotation='vertical',color='white')

#axes.set_title("Landslide Stability fields")
#plt.xscale('log')
plt.xlim(0.25,5)
plt.ylim(0.2,5.8)#plt.ylim(0.1,6) 
#plt.yscale('log')
plt.ylabel(r'frequency ($Hz$)',fontsize=22)
plt.xlabel(r"incident PGA ($m/s^{2}$)",fontsize=22)
plt.yticks(fontsize=19)
plt.xticks(fontsize=19)
plt.text(0.01,0.92,'b)',transform=axes.transAxes,fontsize=22,color='white')
plt.text(0.85,0.1,'40 m',transform=axes.transAxes,fontsize=20,color='white')
plt.text(0.85,0.01,'a < b',transform=axes.transAxes,fontsize=20,color='white')
#plt.xlabel(r"velocity ($m/s$)")

# plt.plot(0.45,1.7,marker='o',mec='white',mfc='indianred',ms=16,mew=2)
# plt.plot(0.5,1.6,marker='o',mec='white',mfc='indianred',ms=16,mew=2)
# plt.plot(1.2,1.6,marker='o',mec='white',mfc='indianred',ms=16,mew=2)
# plt.plot(1.6,1.25,marker='o',mec='white',mfc='royalblue',ms=16,mew=2) #royalblue
# plt.plot(2.4,1.3,marker='o',mec='white',mfc='royalblue',ms=16,mew=2) #royalblue

param={}
param['cs']=par['cs1'][0]
param['rho']=par['rho1'][0] #2000
param['mu']= par['cs1'][0]**2 * par['rho1'][0]#5*1e9
param['h'] = par['h'][0]#-1 #-5
param['tau_c'] = 60200#75000#165424#56200#26311#115424.70790786108 #81000 #45000 #26311 #123118
f=np.linspace(0.25,10,400)
analy=pganalytichugo(f,param)
plt.xlim(0.05,5)
plt.ylim(0.25,6)
plt.plot(pganalytichugo(f,param),f,linestyle='dashed',color='yellow', lw=5)
#plt.show()


plt.savefig('svw_pga40_stability',bbox_inches='tight',dpi=1000,format='pdf')
plt.show()
###








##pganalyfull :
param={}
param['cs']=par['cs1'][0]
param['rho']=par['rho1'][0] #2000
param['mu']= par['cs1'][0]**2 * par['rho1'][0]#5*1e9
param['h'] = par['h'][0]-5#-1 #-5
param['tau_c'] = 115424#56200#26311#115424.70790786108 #81000 #45000 #26311 #123118
param['sig0'] = -781813.5990624019#-781813.5990624019#par['Tn0'][0]
param['a'] = par['A'][0]
param['b'] = par['B'][0]
param['theta0'] = par['theta'][0]
param['v0'] = par['V0'][0]
param['Dc'] = par['Dc'][0]
param['duration']=100
f=np.linspace(0.25,10,400)
#plt.text(0.28,0.84,"1000",fontsize=16,rotation='vertical',color='yellow')
#plt.text(0.5,0.84,"100",fontsize=16,rotation='vertical',color='yellow')
#plt.xlim(0.05,5)
#plt.ylim(0,6.8)
plt.plot(pganalytihugo_full(f,param),f,linestyle='dashed',color='yellow', lw=5)

plt.plot(pganalytihugo_full_thetabis(f,param),f, lw=5)
plt.xlim(0,10000000)
plt.show()


##pga
param['duration']=20
ctepga=param['a']*10/2*(np.log(np.pi*(param['Dc']/(param['v0']*param['duration']))**2)-np.log(np.log(np.pi*(param['Dc']/(param['v0']*param['duration']))**2)))
plt.plot(np.linspace(ctepga,ctepga,2),np.linspace(0,6,2),linestyle='dashed',color='red', lw=5)

#pgv
#param['duration']=2000
f=np.linspace(0.25,10,400)
ctepga=[]
for i in range(len(f)):
    x= f[i] * ( np.pi * 2 )
    z= complex(0,x)
    #ctepga.append(np.pi*2*f[i]*np.abs(param['a']*param['sig0']*param['cs']/((1 - np.exp(-z*2*param['h']/param['cs']))*param['mu']))*(np.log(np.pi*(param['Dc']/(param['v0']*param['duration']))**2)-np.log(np.log(np.pi*(param['Dc']/(param['v0']*param['duration']))**2))))
    ctepga.append(np.pi*2*f[i]*np.abs(param['a']*param['sig0']*param['cs']/(2*param['mu']))*(np.log(np.pi*(param['Dc']/(param['v0']*param['duration']))**2)-np.log(np.log(np.pi*(param['Dc']/(param['v0']*param['duration']))**2))))
plt.plot(ctepga,f,linestyle='dashed',color='cyan', lw=5)

plt.plot(np.linspace(ctepga,ctepga,2),np.linspace(0,6,2),linestyle='dashed',color='cyan', lw=5)

plt.plot(np.linspace(ctepga,ctepga,2),np.linspace(0,6,2),linestyle='dashed',color='yellow', lw=5)



##
param={}
param['cs']=par['cs1'][0]
param['rho']=par['rho1'][0] #2000
param['mu']= par['cs1'][0]**2 * par['rho1'][0]#5*1e9
param['h'] = par['h'][0]-5#-1 
param['tau_c'] = 115424#56200#26311#115424.70790786108 #81000 #45000 #26311 #123118
param['sig0'] = -781813.5990624019#par['Tn0'][0]
param['a'] = par['A'][0]
param['b'] = par['B'][0]
param['theta0'] = par['theta'][0]
param['v0'] = par['V0'][0]
param['Dc'] = par['Dc'][0]


param['duration']=200
f=np.linspace(0.0025,20,800)
brkpnt=param['cs']/(2*param['h']*np.pi)
truc=pganalytihugo_full(f,param)
truc2=[x/(param['a']*10) for x in truc]
fig, axes = plt.subplots(nrows=1, figsize=(7, 5), constrained_layout=True)
cmap=plt.cm.seismic
clrs=cmap(np.linspace(0,1,2))
#red
for i in range(1000):
    plt.plot([2,100],[-0.2+i*0.008,-0.2+i*0.008],color=clrs[1],lw=3) #'red'
#blue
for i in range(1000):
    truc2n=[x+i*0.2 for x in truc2]
    plt.plot(truc2n,f/(param['cs']/(2*param['h'])),color=clrs[0], lw=3) #'blue'
#real
plt.plot(truc2,f/(param['cs']/(2*param['h'])),color='k', lw=3)
plt.plot([2,100],[1,1],linestyle='dashed',alpha=0.66,color='grey',lw=1.5,label=r"$\frac{nC_s}{2 h}/\frac{C_s}{2 h}$")
plt.plot([2,100],[2,2],linestyle='dashed',alpha=0.66,color='grey',lw=1.5)
plt.plot([2,100],[3,3],linestyle='dashed',alpha=0.66,color='grey',lw=1.5)
plt.plot([2,100],[4,4],linestyle='dashed',alpha=0.66,color='grey',lw=1.5)
plt.plot([2,100],[brkpnt/(param['cs']/(2*param['h'])),brkpnt/(param['cs']/(2*param['h']))],linestyle='dotted',alpha=0.66,color='black',lw=1.5,label=r"$\frac{C_s}{2 h \pi}/\frac{C_s}{2 h}$")
plt.xlim(0.25/(param['a']*10),10/(param['a']*10))
plt.ylim(0,4)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.ylabel(r'normalised frequency '+r"$\frac{f 2 h}{C_s}$",fontsize=14)
plt.xlabel(r"normalised incident PGA "+r"$\frac{A_{inc}}{a g}$",fontsize=14)
plt.legend(fontsize=12)


param['duration']=2000
colorpga='darkgreen' #limegreen forestgreen darkgreen
ctepga=param['a']*10/2*(np.log(np.pi*(param['Dc']/(param['v0']*param['duration']))**2)-np.log(np.log(np.pi*(param['Dc']/(param['v0']*param['duration']))**2)))
plt.plot(np.linspace(ctepga/(param['a']*10),ctepga/(param['a']*10),2),np.linspace(0,brkpnt/(param['cs']/(2*param['h'])),2),color=colorpga, lw=3) 

#pgv
#param['duration']=2000
f=np.linspace(brkpnt,20,800)
ctepga=[]
for i in range(len(f)):
    x= f[i] * ( np.pi * 2 )
    z= complex(0,x)
    #ctepga.append(np.pi*2*f[i]*np.abs(param['a']*param['sig0']*param['cs']/((1 - np.exp(-z*2*param['h']/param['cs']))*param['mu']))*(np.log(np.pi*(param['Dc']/(param['v0']*param['duration']))**2)-np.log(np.log(np.pi*(param['Dc']/(param['v0']*param['duration']))**2))))
    ctepga.append(np.pi/(param['a']*10)*2*f[i]*np.abs(param['a']*param['sig0']*param['cs']/(2*param['mu']))*(np.log(np.pi*(param['Dc']/(param['v0']*param['duration']))**2)-np.log(np.log(np.pi*(param['Dc']/(param['v0']*param['duration']))**2))))
plt.plot(ctepga,f/(param['cs']/(2*param['h'])),color=colorpga, lw=3)

plt.savefig('svw_pga40_stability_5',bbox_inches='tight',dpi=300,format='png')
plt.show()

# fig.canvas.draw()
# plt.xlim(0.25,10)
# #plt.ylim(0.053,5.3)#
# plt.ylim(0,4)#plt.ylim(0.1,6) 
# #plt.yscale('log')
# plt.ylabel(r'normalised frequency')
# plt.xlabel(r"incident PGA ($m/s^{2}$)")
# f=np.linspace(0.25,20,800)
# truct=[0,2.5,5,7.5,10,12.5,15,17.5,20]
# flabel=[x/(param['cs']/(2*param['h'])) for x in truct]
# #plt.text(0.28,0.84,"1000",fontsize=16,rotation='vertical',color='yellow')
# #plt.text(0.5,0.84,"100",fontsize=16,rotation='vertical',color='yellow')
# #plt.xlim(0.05,5)
# #plt.ylim(0,6.8)
# truc=pganalytihugo_full(f,param)
# truc2=[x/(param['a']*10) for x in truc]
# truc=[str(x/(param['cs']/(2*param['h']))) for x in truct]
# plt.plot(pganalytihugo_full(f,param),f,color='k', lw=3)
# plt.plot(truc/(param['a']*10),f/(param['cs']/(2*param['h'])),color='k', lw=3)
# labels=[item.get_text() for item in axes.get_yticklabels()]
# labels=['0','0.6','1.3','2','2.6','3.3','4','4.6','5.3']
# axes.set_yticklabels(labels)
# #plt.yticks(flabel)





### subplot résumant expérience
fig, axes = plt.subplots(nrows=2,ncols=2, figsize=(20, 10), constrained_layout=True)
i=15
exptype=1

cache=0
for j in range(5): #len(listd2)-1
    if j<10: cnbr='000'+str(j)
    elif j>=10 and j<100: cnbr='00'+str(j)
    elif j>=100 and j<1000: cnbr='0'+str(j)
    else: cnbr=str(j)
    file_src=dir+cn+'/'+cn_1[0]+'/'+ld[i]+'/'+cnbr
par=pickle.load(open(dir+cn+'/'+cn_1[0]+'/'+ld[i]+'/par.pkl','rb'))
stf=pickle.load(open(file_src+'/stf.pkl','rb'))
data=pickle.load(open(file_src+'/data.pkl','rb'))
#VLSdiff=pickle.load(open(dir+cn+'_resultplots_'+exptype'/vlsdiff.pkl','rb'))

# fig 1. haut gauche - onde input + params
#disp=[0]
#[disp.append(disp[j-1]+stf['stf'][j+1]*(stf['t'][j+1]-stf['t'][j])) for j in range(1,len(stf['t'])-2) ]
axes[0,0].plot(stf['t'],stf['stf'],'red')
axes[0,0].set_ylabel("Velocity (m/s)")
#ax2=axes[0].twinx()
# make a plot with different y-axis using second axis object
#ax2.plot(stf['t'][:-1],disp,color="blue")
#ax2.set_ylabel("Displacement (m)",color="blue",fontsize=14)
axes[0,0].set_title("Input wave")
textstr = '\n'.join((
    r'$\mathrm{ampli}_{\mathrm{max}}=%.2e$ m/s ' % (par['ampli'][j], ),
    r'$\mathrm{f}=%.2e$ Hz' % (par['f0'][i], ), 
    r'$\sigma=%.2e$' % (par['std'][i], ),
    r'$\mathrm{H}=%.2e$ m' % (par['h'][i], ),
    r'$\mathrm{A - B}=%.2e$' % (par['A'][i]-par['B'][i], ),
    r'$\mathrm{State Law}= Slip Law',
    r'$\mathrm{Angle}=%.2e$ rad' % (par['angle'][i], ),
    r'$\rho=%.2e$ kg/m3' % (par['rho1'][i], ),
    r'$V_{0}=%.2e$ m/s' % (par['V0'][i], )))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
axes[0,0].text(0.05, 0.95, textstr, transform=axes[0,0].transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

# fig 2. haut droite - onde output 
Vnorm=data['v']/par['V0'][i]
timex=np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache)
axes[0,1].plot(np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache),Vnorm[cache:,int(data['nx']/2)],'blue')
axes[0,1].set_xlabel("time (s)")
axes[0,1].set_ylabel("Velocity normalized", color='blue', fontsize=14)
axes[0,1].set_yscale('log')
ax2=axes[0,1].twinx()
ax2.plot(np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache),data['d'][cache:,int(data['nx']/2)],'grey',lw=4)
ax2.set_ylabel('Displacement (m)', color='grey',fontsize=14)
axes[0,1].set_title("Velocity middle fault -- ")
textstr = '\n'.join((
    r'$ \mathrm{Last} V / V_{0} =$%.2e ' % (np.mean(data['v'][-10000:-1,int(data['nx']/2)])/par['V0'][i], ),
    r'Acceleration $=%.2e$' % ((data['v'][-1,int(data['nx']/2)]-data['v'][-10000,int(data['nx']/2)])/(timex[-1]-timex[-10000])) ))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
axes[0,1].text(0.05, 0.95, textstr, transform=axes[0,1].transAxes, fontsize=14,
        verticalali gnment='top', bbox=props)

# fig 3. bas gauche - tau output 
axes[1,0].plot(np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache),data['st'][cache:,int(data['nx']/2)]/par['Tt0'][i],'orange')
axes[1,0].set_xlabel("time (s)")
axes[1,0].set_ylabel("Tangential stress normalized", color='orange', fontsize=14)
axes[1,0].set_title("Tangential stress")

# fig 4. bas droit - VLS diff
axes[1,1].plot(np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache),((data['st'][cache:,int(data['nx']/2)])/par['Tn0'][i]+np.log(Vnorm[cache:,int(data['nx']/2)])*par['A'][i]),'black',lw=2)
axes[1,1].set_xlabel("time (s)")
axes[1,1].set_ylabel("Écart", color='black', fontsize=14)
axes[1,1].set_title("Écart à la solution analytique")

plt.savefig('exp-res_'+cn+'_'+cn_1[0]+str(i)+'_'+str(cnbr),bbox_inches='tight',dpi=300,format='png')
plt.show()




############ comparison sin vs gauss*sin

def gaussxsin(ampli,f0,std,onset,t):
    #sinarg=[math.radians(f0*np.pi*2*(x-onset)) for x in t]
    sinarg=f0*np.pi*2*(t-onset)
    return ampli*np.exp(-1/2*((t-onset)/std)**2)*np.sin(sinarg)

def sin(ampli,f0,onset,t):
    sinarg=f0*np.pi*2*(t-onset)
    return ampli*np.sin(sinarg)

f0=0.5
ampli=0.3*10/(2*np.pi*f0)#par['ampli'][44]
std=600
onset=30
t=np.linspace(0,60,379063)

gauss=gaussxsin(ampli,f0,std,onset,t)
sinc=sin(ampli,f0,onset,t)

stftot=gauss
hcs=45*2/300
superpo1 = list(np.zeros(np.int(hcs/(t[1]-t[0])+1)))
superpo1.extend(list(stftot))
superpo2 = list(stftot)
superpo2.extend(list(np.zeros(np.int(hcs/(t[1]-t[0])+1))))
superpototg=[]
[superpototg.append(superpo2[i]-superpo1[i]) for i in range(len(superpo2))]
pga_inc_gauss=(2*np.pi*f0)*np.max(superpototg)
tau_inc_gauss=2000*300*np.max(superpototg)

stftot=sinc
superpo1 = list(np.zeros(np.int(hcs/(t[1]-t[0])+1)))
superpo1.extend(list(stftot))
superpo2 = list(stftot)
superpo2.extend(list(np.zeros(np.int(hcs/(t[1]-t[0])+1))))
superpotots=[]
[superpotots.append(superpo2[i]-superpo1[i]) for i in range(len(superpo2))]
pga_inc_sin=(2*np.pi*f0)*np.max(superpotots)
tau_inc_sin=2000*300*np.max(superpotots)

print(tau_inc_gauss)
print(tau_inc_sin)



plt.plot(t,sinc,'k')
plt.plot(t,gauss,'r')
plt.xlim(20,40)
plt.savefig('gauss-v-sin',bbox_inches='tight',dpi=300,format='png')
plt.show()

ttg=np.linspace(0,60,len(superpototg))
tttot=np.linspace(0,60,len(superpotot))

plt.plot(stftot['t'],stftot['stf'],'b')
plt.plot(t,gauss,'r')
plt.xlim(20,40)
plt.savefig('gauss-v-gauss',bbox_inches='tight',dpi=300,format='png')
plt.show()