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
import gc
from function_read import sem2d_read_fault, sem2d_read_STF

dir='' #working directory
cn='SF1_svw_1h15-s-5--35-aw-5' #codename of the exp
cn_1=['anglewave']#'B'#,'f0','std'] #subcodename of the parameter alterated
subexp='angle' #pga
linear_sol=True
cache=0

# verify if dir of output exist
ref_dir=dir+cn+'/'+cn_1[0]
if path.exists(ref_dir):
    print('going in : '+ref_dir)
else:
    print('No directory associated')

exiti=0
listd1=os.listdir(ref_dir)
for i in range(len(listd1)):
    if exiti==1: #2
        break
    listd2=os.listdir(ref_dir+'/'+listd1[i])
    par=pickle.load(open(ref_dir+'/'+listd1[i]+'/par.pkl','rb'))
    nameoutput='/plots/plots_'+cn_1[0]+'='+str(par[cn_1[0]][0])
    if path.exists(dir+cn+nameoutput):
        print(dir+cn+nameoutput+' already exist')
        VLSdiff=pickle.load(open(dir+cn+nameoutput+'/vlsdiff.pkl','rb'))
        whencheck=pickle.load(open(dir+cn+nameoutput+'/whencheck.pkl','rb'))
        norml2tot=pickle.load(open(dir+cn+nameoutput+'/norml2.pkl','rb'))
    else:
        exiti+=1
        if not path.exists(dir+cn+'/plots'):
            sp.run(['mkdir',dir+cn+'/plots'])
        sp.run(['mkdir',dir+cn+nameoutput])
        sp.run(['mkdir',dir+cn+nameoutput+'/velocity'])
        sp.run(['mkdir',dir+cn+nameoutput+'/tengential_stress'])
        sp.run(['mkdir',dir+cn+nameoutput+'/vs_linear-sol'])
        sp.run(['mkdir',dir+cn+nameoutput+'/theta'])
        VLSdiff=[]
        whencheck=[]
        norml2tot=[]
    for j in range(len(listd2)-1):
        if j<10: cnbr='000'+str(j)
        elif j>=10 and j<100: cnbr='00'+str(j)
        elif j>=100 and j<1000: cnbr='0'+str(j)
        else: cnbr=str(j)
        file_src=ref_dir+'/'+listd1[i]+'/'+cnbr
        if path.exists(dir+cn+nameoutput+'/velocity/vel_'+cnbr+'.png'):
            continue
        else:
            try:
                listd3=os.listdir(file_src)
                print('going in : '+file_src)
                if any(e in 'data.pkl' for e in listd3): # check output already written
                    data=pickle.load(open(file_src+'/data.pkl','rb'))
                    source_time_fct=pickle.load(open(file_src+'/stf.pkl','rb'))
                else:
                    print('Creating output in : '+file_src)
                    for k in range(len(listd3)):
                        if listd3[k][0:3]=='Flt':
                            Fltname=listd3[k][0:5]
                    data=sem2d_read_fault(file_src,Fltname)
                    source_time_fct=sem2d_read_STF(file_src)
                    pickle.dump(data,open(file_src+'/data.pkl','wb'),protocol=4)
                    pickle.dump(source_time_fct,open(file_src+'/stf.pkl','wb'),protocol=4)
                    sp.run('rm -r '+file_src+'/*sem2d*',shell=True)
                try:
                    if len(data['v'])>2:
                        Vnorm=data['v']/par['V0'][j]
                        stress_eff=(data['st'][:,int(data['nx']/2)]/par['Tt0'][j]+data['sn'][:,int(data['nx']/2)]/par['Tn0'][j])/2
                        ## fig vitesse norm
                        plt.figure()
                        plt.plot(np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache),Vnorm[cache:,int(data['nx']/2)])
                        plt.xlabel("time (s)")
                        plt.ylabel("Velocity normalized")
                        plt.yscale('log')
                        plt.title("Velocity middle fault -- "+subexp+"="+str("{:.2e}".format(par[subexp][j])))
                        plt.savefig(dir+cn+nameoutput+'/velocity/vel_'+cnbr+'.png',bbox_inches='tight',dpi=300,format='png')
                        plt.close()
                        ## fig stress tangential
                        plt.figure()
                        plt.plot(np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache),data['st'][cache:,int(data['nx']/2)]/par['Tt0'][j])
                        plt.xlabel("time (s)")
                        plt.ylabel("Tangential stress normalized")
                        plt.title("Tangential stress middle fault -- "+subexp+"="+str("{:.2e}".format(par[subexp][j])))
                        plt.savefig(dir+cn+nameoutput+'/tengential_stress/st_'+cnbr+'.png',bbox_inches='tight',dpi=300,format='png')
                        plt.close()
                        ## fig theta
                        plt.figure()
                        plt.plot(np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache),data['theta'][cache:,int(data['nx']/2)]/par['theta'][j])
                        plt.xlabel("time (s)")
                        plt.ylabel("Theta/theta0")
                        plt.title("Theta normalised -- "+subexp+"="+str("{:.2e}".format(par[subexp][j])))
                        plt.savefig(dir+cn+nameoutput+'/theta/theta_'+cnbr+'.png',bbox_inches='tight',dpi=300,format='png')
                        plt.close()
                        if linear_sol:
                            Vnew=Vnorm[cache:,int(data['nx']/2)]
                            mask=np.where(Vnew<0)[0]
                            for m in range(len(mask)):
                                Vnew[mask[m]]=par['V0'][j]*0.1
                            plt.figure()
                            plt.plot(np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache),(-np.log(Vnew)*par['A'][j]),lw=2,label='-log(V/V0)*A')
                            plt.plot(np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache),((data['st'][cache:,int(data['nx']/2)])/par['Tn0'][j]),linestyle='--',lw=1,label='tau/tau0')
                            plt.legend()
                            plt.xlabel("time (s)")
                            plt.ylabel("Tstress normalized ou -Alog(V/V0)")
                            plt.title("Tsress vs. -Alog(V/V0) -- "+subexp+"="+str("{:.2e}".format(par[subexp][j])))
                            plt.savefig(dir+cn+nameoutput+'/vs_linear-sol/VLS_'+cnbr+'.png',bbox_inches='tight',dpi=300,format='png')
                            plt.close()

                            plt.figure()
                            plt.plot(np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache),((data['st'][cache:,int(data['nx']/2)])/par['Tn0'][j]+np.log(Vnew)*par['A'][j]),lw=2)
                            plt.xlabel("time (s)")
                            plt.ylabel("-Alog(V/V0) - Tstress normalized ")
                            plt.title("VLS diff-- "+subexp+"="+str("{:.2e}".format(par[subexp][j])))
                            plt.savefig(dir+cn+nameoutput+'/vs_linear-sol/VLSdiff_'+cnbr+'.png',bbox_inches='tight',dpi=300,format='png')
                            plt.close()


                            ##
                            VLSdiffmean=np.mean(np.mean((data['st'][cache:,int(data['nx']/2)])/par['Tn0'][j]+np.log(Vnew)*par['A'][j]))
                            VLSdiff.append(VLSdiffmean)
                            pickle.dump(VLSdiff,open(dir+cn+nameoutput+'/vlsdiff.pkl','wb'),protocol=4)
                            #whencheck.append(j)
                            pickle.dump(whencheck,open(dir+cn+nameoutput+'/whencheck.pkl','wb'),protocol=4)
                            norml2=np.sqrt(np.sum(((data['st'][cache:,int(data['nx']/2)])/par['Tn0'][j]+np.log(Vnew)*par['A'][j])**2))#/len(Vnew))
                            norml2tot.append(norml2)
                            #testl2=np.linalg.norm((data['st'][cache:,int(data['nx']/2)])/par['Tn0'][j]+np.log(Vnew)*par['A'][j],2)
                            pickle.dump(norml2tot,open(dir+cn+nameoutput+'/norml2.pkl','wb'),protocol=4)
                            ##

                            # plt.figure()
                            # plt.plot(np.linspace(0,data['dt']*(len(data['v'])-cache),len(data['v'])-cache),(-np.log(Vnew)*par['A'][j])/((data['st'][cache:,int(data['nx']/2)])/par['Tn0'][j]),lw=2,label='-log(V/V0)*A / (tau/tau0)')
                            # plt.yscale('log')
                            # plt.legend()
                            # plt.xlabel("time (s)")
                            # plt.ylabel("-Alog(V/V0) / (tau/tau0)")
                            # plt.title("Tsress vs. -Alog(V/V0) -- "+cn_1[0]+"="+str("{:.2e}".format(par[cn_1[0]][j])))
                            # plt.savefig(dir+cn+'_'+nameoutput+'/vs_linear-sol/VLSratio_'+cnbr+'.png',bbox_inches='tight',dpi=300,format='png')
                            # plt.close()
                    else:
                        continue
                except(NameError):
                    continue
                whencheck.append(j)
            except(ValueError):
                continue
    if path.exists(dir+cn+nameoutput+'/velocity/vel_'+cnbr+'.png'):
        continue
    else:
        del listd2
        del par
        del nameoutput
        del VLSdiff
        del norml2tot
        del whencheck
        del data
        del source_time_fct
        del Vnorm
        del stress_eff
        del mask
        del Vnew
        del VLSdiffmean
        del norml2
        gc.collect()

    par=pickle.load(open(ref_dir+'/'+listd1[i]+'/par.pkl','rb'))
    nameoutput='/plots/plots_'+cn_1[0]+'='+str(par[cn_1[0]][0])
    VLSdiff=pickle.load(open(dir+cn+nameoutput+'/vlsdiff.pkl','rb'))
    if linear_sol:
        plt.figure()
        #plt.plot(par[cn_1[0]],VLSratio)
        #del(par[cn_1[0]][whencheck])
        plt.plot(par[subexp][:len(VLSdiff)],VLSdiff)
        plt.xlabel("facteur étudié : "+subexp)
        plt.ylabel("Ratio -Alog(V/V0) / t/t0")
        plt.xscale('log')
        #plt.yscale('log')
        #plt.title("Tsress vs. -Alog(V/V0) -- "+cn_1[0]+"="+str("{:.2e}".format(par[cn_1[0]][j])))
        plt.savefig(dir+cn+nameoutput+'/VLSratioall-fac-'+cn+'.png',bbox_inches='tight',dpi=300,format='png')
        plt.close()

        plt.figure()
        #plt.plot(VLSratio)
        plt.plot(VLSdiff)
        plt.xlabel("Exp number")
        plt.ylabel("Ratio -Alog(V/V0) / t/t0")
        plt.savefig(dir+cn+nameoutput+'/VLSratioall-exp-'+cn+'.png',bbox_inches='tight',dpi=300,format='png')
        plt.close()
        pickle.dump(VLSdiff,open(dir+cn+nameoutput+'/vlsdiff.pkl','wb'),protocol=4)