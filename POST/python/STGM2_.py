"""
Created on Thu Jan 13 14:38:33 2022

Programme STGM2-plot.py
    FUNCTION models_path=get_path_STGM2()
        !!!toujours active!!!
        effectue via clavier la recherche du modele numerique
        retourne le path de toutes les simulations pour ce modele
        
    FUNCTION all_(models_path,ref):
        plot toutes les figure pour la variable : ref
        :ref=['contour','Run_speed','Run_acc','G0','Peak_sliprate','slip','stress_drop','sliprate','stress_change']
    
    FUNCTION all__(models_path,ref1,ref2):
        plot la superposition de deux variables : ref1,ref2
        if
        :ref1 and ref2 in [Run_speed','Run_acc','G0','Peak_sliprate','slip','stress_drop']
        :ref1 and ref2 in['sliprate','stress_change']
        
        
        

@author: florian llorens
"""

# -*- coding: utf-8 -*-
#param : [m]
#    '-1'                   le dernier de la liste
#    'm=[0 1 2 3 ...]'      le m
#
#data :
#nx      = data['nx']                #nx
#nt      = data['nt']                #nt  
#dt      = data['dt']                #dt
#x       = data['x']/1e3        #km
#v       = data['v']            #sliprate
#d       = data['d']            #d ?
#st      = data['st']           #st stress change
#sn      = data['sn']           #sn stress change ?
#mu      = data['mu']           #coeff de friction ? de lame ?
#st0      = data['st0']         #st0
#sn0      = data['sn0']         #sn0
#mu0      = data['mu0']         #mu0  
#
#x=[]
#sliprate=[]
#stress_change=[]
#
#xx=[]
#Heal_t=[]
#Rup_speed=[]
#Rup_acc=[]
#Gc=[]
#G0=[]
#Peak_sliprate=[]
#slip=[]
#stress_drop=[]





import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from os import listdir,mkdir
from os.path import isfile,isdir,join,exists


#------------------------------------------------------------------------------------
#                                '''''''path'''''''
#------------------------------------------------------------------------------------
local_path='D:/STAGEM2/output/'#'/u/moana/user/llorens/STGM2/output/'#dossier in  
out_path='D:/STAGEM2/scripts/'#"/u/moana/user/llorens/STGM2/scripts/"#dossier out
local_process_path='D:/STAGEM2/scripts/'#"/u/moana/user/llorens/STGM2/scripts/"#dossiers de la source
alpha=1/5#Trigger
#------------------------------------------------------------------------------------
#                                '''''''saisie clavier du repertoire'''''''
#------------------------------------------------------------------------------------

def get_path_STGM2():    
    global local_path#="D:/STAGEM2/output/"#dossiers des modeles de simulation
    m=[]
    ch=-1#choix du fichier modele
    fichiers = [f for f in listdir(local_path) if not isfile(join(local_path, f))]
    def ord_path(models):
        models_code=[];ord_mod=[]
        for m in models:
            code=''
            for c in m:
                if c!='_':code+=c
                else:break
            models_code.append(int(code))
        models_code=sorted(models_code)
        for N in models_code:
            for m in models:
                code_=''
                for c in m:
                    if c!='_':code_+=c
                    else:break
                if N==int(code_):ord_mod.append(m)
        return ord_mod
    i=0
    for fi in fichiers:
        i=i+1
        m.append(fi);print('\n',i,'>',fi)
    while ch not in list(range(1,1+len(fichiers))):
        ch=input('>')
        try:   ch=int(ch)
        except ValueError:print('!error : function get_path_STGM2()')
    _path=local_path+m[ch-1]
    models=[f for f in listdir(_path) if not isfile(join(_path, f))]#toutes les simu d un modele
    models=ord_path(models)
    i=0
    for mod in models:
        models[i]=_path+'/'+mod
        #print('>',mod)
        i=i+1
    return models#contient les path de chaque simu disponible




#------------------------------------------------------------------------------------
#                                ''''''PROCESS''''''''''
#------------------------------------------------------------------------------------

def sem2d_read_fault(model_name,fault_name):
    
    # length of the tag at the begining and end of a binary record
    # in number of single precision words (4*bytes)
    LENTAG = 2; # gfortran older versions
    LENTAG = 1;
    
    # assumes header file name is FltXX_sem2d.hdr
    if not isdir(model_name):
        print("Wrong path to the model directory...")
        quit()
    headfile_exist = isfile(model_name+"/"+fault_name+"_sem2d.hdr")
    initfile_exist = isfile(model_name+"/"+fault_name+"_init_sem2d.tab")
    datafile_exist = isfile(model_name+"/"+fault_name+"_sem2d.dat")
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

#------------------------------------------------------------------------------------
#                                ''''''PLOT''''''''''
#------------------------------------------------------------------------------------


def get_snapshot_data(path,nb):
    snap=open(path+"/snapshot_"+nb+".dat","r")
    fic=snap.readlines()    
    X=[];Sliprate=[];Stress_change=[]
    for li in fic:
        d=li.split('  ')
        X.append(float(d[0]))
        Sliprate.append(float(d[1]))
        Stress_change.append(float(d[2]))
    snap.close()
    return X,Sliprate,Stress_change

def get_asv_data(path):
    al_s=open(path+"/along_strike_values.dat","r")
    fic=al_s.readlines()  
    X=[];Rup_t=[];Heal_t=[];Rup_speed=[];Rup_acc=[];Gc=[];G0=[];Peak_sliprate=[];slip=[];stress_drop=[]
    for li in fic:
        d=li.split('  ')
        X.append(float(d[0]))
        Rup_t.append(float(d[1]))
        Heal_t.append(float(d[2]))
        Rup_speed.append(float(d[3]))
        Rup_acc.append(float(d[4]))
        Gc.append(float(d[5]))
        G0.append(float(d[6]))
        Peak_sliprate.append(float(d[7]))
        slip.append(float(d[8]))
        stress_drop.append(float(d[9][:-1]))
    al_s.close()
    return X,Rup_t,Heal_t,Rup_speed,Rup_acc,Gc,G0,Peak_sliprate,slip,stress_drop
 
def get_max_contour(data):
    X= data['x']/1e3;V= data['v']
    NX=data['nx'];NT=data['nt'];DT=data['dt']
    global alpha
    x , t = np.meshgrid(np.linspace(0,X[-1],NX),np.linspace(0,NT*DT,NT))
    x_max=-1;t_max=-1;i=0;threshold=-1;v_avg=[]
    for li in np.transpose(V):#sur l espace
        v_avg_=sum(li)/len(li)
        v_avg.append(v_avg_)
        if v_avg_>=threshold:threshold=v_avg_
        elif v_avg_>(threshold*alpha/10) and v_avg_< (threshold*alpha):x_max=X[i]
        i=i+1
    i=0
    for li in V:#sur le temps
        v_avg_=sum(li)/len(li)
        if v_avg_>=threshold:threshold=v_avg_
        elif v_avg_>(threshold*alpha/10) and v_avg_< (threshold*alpha):t_max=DT*i+DT
        i=i+1
    return x_max,t_max

def get_id(liste,rf):    
    i=0
    for elem in liste: 
        if str(elem)==str(rf):return i
        i=i+1
    return -1

#-----------------------------------------------------------------------------------------

def one_anim(data,fr,ref):
    def av_anim(j):
        Y=f[j,:];line.set_data(X,Y);return line,
    X = data['x']/1e3;v = data['v'];st = data['st'];sn = data['sn'];mu=data['mu'];d = data['d']   
    fig= plt.figure()
    ax=plt.axes();ax.set_xlim(0,np.max(X));ax.set_xlabel('X')
    line, = ax.plot([], []) 
    f=[]
    if ref=='v':f=v;ax.set_ylim(0,np.max(v));ax.set_ylabel('v')    
    elif ref=='st':f=st;ax.set_ylim(0,np.max(st));ax.set_ylabel('st')
    elif ref=='sn':f=sn;ax.set_ylim(0,np.max(sn));ax.set_ylabel('sn')
    elif ref=='mu':f=mu;ax.set_ylim(0,np.max(mu));ax.set_ylabel('mu')
    elif ref=='d':f=d;ax.set_ylim(0,np.max(d));ax.set_ylabel('d')
    ax.plot(X,f[0,:],'r') 
    animation.FuncAnimation(fig, av_anim, frames=fr, blit=True, interval=20, repeat=False)

def plot_1(X,y,ref,titre):
    fig, ax1 = plt.subplots()
    color = 'tab:orange'
    plt.title(titre)
    ax1.set_xlabel('X (Km)');ax1.set_ylabel(ref, color=color)
    ax1.plot(X,y, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    plt.show()

def plot_2(X,y1,y2,ref1,ref2,titre):
    fig, ax1 = plt.subplots()
    color = 'tab:orange'
    plt.title(titre)
    ax1.set_xlabel('X (Km)');ax1.set_ylabel(ref1, color=color)
    ax1.plot(X,y1, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:blue'
    ax2.set_ylabel(ref2, color=color)  # we already handled the x-label with ax1
    ax2.plot(X,y2, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()
    
def plot_2_xlab(X,y1,y2,ref1,ref2,titre,xlab):
    fig, ax1 = plt.subplots()
    color = 'tab:orange'
    plt.title(titre)
    ax1.set_xlabel(xlab);ax1.set_ylabel(ref1, color=color)
    ax1.plot(X,y1, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:blue'
    ax2.set_ylabel(ref2, color=color)  # we already handled the x-label with ax1
    ax2.plot(X,y2, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()
    
def plot_contour_v(data,title):
    X= data['x']/1e3;V= data['v'];NX=data['nx'];NT=data['nt'];DT=data['dt']
    x , t = np.meshgrid(np.linspace(0,X[-1],NX),np.linspace(0,NT*DT,NT))     
    fig=plt.figure()
    ax=fig.add_subplot(111)
    tit=''
    titl=reversed(title)
    for e in titl:
        if e!='/':tit+=e
        else:break
    plt.title("".join(reversed(tit)))
    plt.contourf(x,t,V,25)
    plt.xlabel('X');plt.ylabel('t')   
    txt=r'$y = e^{x}$'#LATEX
    plt.text(0.5,1.09,txt,horizontalalignment='center',verticalalignment='center', transform = ax.transAxes,color='green')
    plt.colorbar()
 
def plot_max_contour(data):
    X= data['x']/1e3;V= data['v'];DT=data['dt']
    global alpha
    x_max=[];v_avg=[];t_max=0;i=0;threshold=-1
    for li in np.transpose(V):#sur l espace
        v_avg_=sum(li)/len(li)
        v_avg.append(v_avg_)
        if v_avg_>=threshold:threshold=v_avg_
        elif v_avg_>(threshold*alpha/10) and v_avg_< (threshold*alpha):x_max=X[i]
        i=i+1
    i=0
    for li in V:#sur le temps
        v_avg_=sum(li)/len(li)
        if v_avg_>=threshold:threshold=v_avg_
        elif v_avg_>(threshold*alpha/10) and v_avg_< (threshold*alpha):t_max=DT*i+DT
        i=i+1
    fig, ax1 = plt.subplots()
    plt.title('mqx contour')
    ax1.set_xlabel('X (Km)');ax1.set_ylabel('avg<Sliprate> ');plt.xlim([0,max(X)])
    plt.title('tmax='+str(t_max)+' xmax :'+str(x_max))
    ax1.plot(X,v_avg)
    ax1.tick_params(axis='y')
    plt.show()

def plot_xmaxtmax_of_dBarr(models_path):
    dBarr=[];x_fin=[];t_fin=[]
    global alpha
    for mod in models_path:
       D=get_Par_data(mod)
       D=convertir_Par_data(D)
       data=sem2d_read_fault(mod,fault_name="Flt01")
       d=D['theta'][1]-D['theta'][0]
       xf,tf=get_max_contour(data) 
       dBarr.append(d);x_fin.append(xf);t_fin.append(tf)
    plot_2(dBarr,x_fin,t_fin,'xf','tf','distance d arret = f(dBarr)')
    
    
def plot_xmaxtmax_of_DcBarr(models_path):
    li_=[];x_fin=[];t_fin=[];penetration=[]
    global alpha
    for mod in models_path:
       D=get_Par_data(mod)
       D=convertir_Par_data(D)
       data=sem2d_read_fault(mod,fault_name="Flt01")
       xf,tf=get_max_contour(data)
       penetration.append(xf-(D['Tt'][1])/1000)
       if type(D['Dc'])==float:li_.append(D['Dc'])
       else :li_.append(D['Dc'][-1]);
       x_fin.append(xf);t_fin.append(tf)
    #plot_2_xlab(li_,x_fin,t_fin,'xf','tf','distance d arret = f(thetaBarr)','Dc')
    if type(D['Dc'])==float:plot_2_xlab(li_,penetration,t_fin,'penetration','tf','penetration barriere= f(*Dc*)','Dc')
    else:plot_2_xlab(li_,penetration,t_fin,'penetration','tf','penetration barriere= f(DcBarr)','Dc')

    
def all_(models_path,ref):
    for m in models_path:
        if ref=='contour':data = sem2d_read_fault(m,fault_name="Flt01");plot_contour_v(data,title=m)
        elif ref=='max_contour':data = sem2d_read_fault(m,fault_name="Flt01");plot_max_contour(data)
        else:
            xx,Rup_t,Heal_t,Rup_speed,Rup_acc,Gc,G0,Peak_sliprate,slip,stress_drop=get_asv_data(m)
            x,sliprate,stress_change=get_snapshot_data(m,nb='3')
            if ref=='Rup_speed':plot_1(xx,Rup_speed,ref,m);
            elif ref=='Rup_acc':plot_1(xx,Rup_acc,ref,m)
            elif ref=='Gc':plot_1(xx,Gc,ref,m)
            elif ref=='G0':plot_1(xx,G0,ref,m)
            elif ref=='Peak_sliprate':plot_1(xx,Peak_sliprate,ref,m)
            elif ref=='slip':plot_1(xx,slip,ref,m)
            elif ref=='stress_drop':plot_1(xx,stress_drop,ref,m)
            elif ref=='sliprate':plot_1(x,sliprate,ref,m)
            elif ref=='stress_change':plot_1(x,stress_change,ref)
        
           
def all__(models_path,ref1,ref2):
    tx_x=['x','sliprate','stress_change']
    tx_xx=['xx','Rup_t','Heal_','Rup_speed','Rup_acc','Gc','G0','Peak_sliprate','slip','stress_drop']
    if ref1 in tx_x and ref2 in tx_x:
        for m in models_path:
            x,sliprate,stress_change=get_snapshot_data(m,nb='3')
            li_x=[ x,sliprate,stress_change]  
            i1=get_id(tx_x,ref1)
            i2=get_id(tx_x,ref2)
            #print(i1,i2,m)
            y1=li_x[i1]
            y2=y2=li_x[i2]
            X=li_x[0]
            plot_2(X,y1,y2,ref1,ref2,titre=m)            
    elif ref1 in tx_xx and ref2 in tx_xx:
        for m in models_path:
            xx,Rup_t,Heal_t,Rup_speed,Rup_acc,Gc,G0,Peak_sliprate,slip,stress_drop=get_asv_data(m)
            li_xx=[xx,Rup_t,Heal_t,Rup_speed,Rup_acc,Gc,G0,Peak_sliprate,slip,stress_drop] 
            i1=get_id(tx_xx,ref1)
            i2=get_id(tx_xx,ref2)
            #print(i1,i2,m)
            y1=li_xx[i1]
            y2=y2=li_xx[i2]
            X=li_xx[0]
            plot_2(X,y1,y2,ref1,ref2,titre=m)
    elif ref1=='contour' or ref2=='contour':
        for m in models_path:
            data = sem2d_read_fault(m,fault_name="Flt01")
            plot_contour_v(data,title=m)
            plot_max_contour(data)




#------------------------------------------------------------------------------------
#                                ''''''PAR.INP'''''''''
#------------------------------------------------------------------------------------
#-----------------------RSF-------------------------------------
# ! BC_DYNFLT_RSF: rate and state dependent friction for dynamic faults
    # ! ARG: Dc       [dble] [0.5d0] Critical slip 
    # ! ARG: MuS      [dble] [0.6d0] Static friction coefficient
    # ! ARG: a        [dble] [0.01d0] Direct effect coefficient
    # ! ARG: b        [dble] [0.02d0] Evolution effect coefficient
    # ! ARG: Vstar    [dble] [1d0] Characteristic or reference slip velocity
    # ! ARG: theta    [dble] [1d0] State variable
    # ! ARG: kind     [int] [1] Type of slip weakening function:
        # !                       1 = #Strong velocity-weakening at high speed  as in Ampuero and Ben-Zion 
        # !                       2 = #Kaneko et al. (2008) Eq. 15 (regularize at v=0 as per Laupsta et al. (2000))
        # !                       3 = #Kaneko et al. (2008) Eq. 15 (regularize at v=0 as per Laupsta et al. (2000))r
        # !                       4 = #????????????????
#-----------------------SWF-------------------------------------
    # ! BC_DYNFLT_SWF: slip weakening friction for dynamic faults
    # ! ARG: kind     [int] [1] Type of slip weakening function:
    # !                       1 = linear
    # !                       2 = exponential
    # ! ARG: Dc       [dble] [0.5d0] Critical slip 
    # ! ARG: MuS      [dble] [0.6d0] Static friction coefficient
    # ! ARG: MuD      [dble] [0.5d0] Dynamic friction coefficient
    # ! ARG: alpha    [dble] [0.0d0] Roughness drag coefficient (normalized by Tn)
    # ! ARG: healing  [log] [F] Instantaneous healing upon slip arrest
    # !               Healing is currently valid only with the leapfrog time scheme
    # !
#---------


def print__(D):
    print('GENERAL       ','ngll     = ',D['ngll'])
    print('GENERAL       ','fmax     = ',D['fmax'])
    print('GENERAL       ','W        = ',D['W'])
    print('GENERAL       ','ndof     = ',D['ndof'])
    print('MESH_CART   ',D['xlim'],D['zlim'],D['nelem'])
    print('MAT_ELASTIC   ','rho      = ',D['rho'])
    print('MAT_ELASTIC   ','cp       = ',D['cp'])
    print('MAT_ELASTIC   ','cs       = ',D['cs'])
    if D['friction']=='RSF':
        print('BC_DYNFLT_RSF ','kind     = ',D['kind'])
        print('BC_DYNFLT_RSF ','Dc       = ',D['Dc'])
        print('BC_DYNFLT_RSF ','MuS      = ',D['MuS'])
        print('BC_DYNFLT_RSF ','a        = ',D['a'])
        print('BC_DYNFLT_RSF ','b        = ',D['b'])
        print('BC_DYNFLT_RSF ','Vstar    = ',D['Vstar'])
        print('BC_DYNFLT_RSF ','Vc       = ',D['Vc'])
        print('BC_DYNFLT_RSF ','theta    = ',D['theta'])
    if D['friction']=='SWF':
        print('BC_DYNFLT_SWF ','kind     = ',D['kind'])
        print('BC_DYNFLT_SWF ','Dc       = ',D['Dc'])
        print('BC_DYNFLT_SWF ','MuS      = ',D['MuS'])
        print('BC_DYNFLT_SWF ','MuD        = ',D['a'])
    print('BC_DYNFLT     ','friction = ',D['friction'])
    print('BC_DYNFLT     ','Tn       = ',D['Tn'])
    print('BC_DYNFLT     ','Tt       = ',D['Tt'])
    print('BC_DYNFLT     ','V        = ',D['V'])

def print_(D):
    print('ngll     = ',D['ngll'])
    print('fmax     = ',D['fmax'])
    print('W        = ',D['W'])
    print('ndof     = ',D['ndof'])
    print('MESH_CART=',D['xlim'],D['zlim'],D['nelem'])
    print('rho      = ',D['rho'])
    print('cp       = ',D['cp'])
    print('cs       = ',D['cs'])
    if D['friction']=='RSF':
        print('kind     = ',D['kind'])
        print('Dc       = ',D['Dc'])
        print('MuS      = ',D['MuS'])
        print('a        = ',D['a'])
        print('b        = ',D['b'])
        print('Vstar    = ',D['Vstar'])
        print('Vc       = ',D['Vc'])
        print('theta    = ',D['theta'])
    if D['friction']=='SWF':
        print('kind     = ',D['kind'])
        print('Dc       = ',D['Dc'])
        print('MuS      = ',D['MuS'])
        print('MuD      = ',D['a'])
    print('friction = ',D['friction'])
    print('Tn       = ',D['Tn'])
    print('Tt       = ',D['Tt'])
    print('V        = ',D['V'])
 
def get_value(chaine,ref):
    #print('get_v:',chaine,ref)
    i=-1
    i=str(chaine).find(ref+'=')
    if i==-1:
        i=str(chaine).find(ref+'H=')
    chaine_=''
    for c in chaine[i:]:
        if c!='/' and c!=',' and c!='\n':
            chaine_=chaine_+c
        else:
            break

    return chaine_[len(ref)+1:]
def update_texte(chaine,ref,value):
    i=-1
    i=str(chaine).find(ref+'=')+len(ref)+1
    if i==-1:
        i=str(chaine).find(ref+'H=')+len(ref)+2
    j=0
    for c in chaine[i:]:
        if c==',' or c=='/' or c=='\n':
            lon=j
            break
        j=j+1
    chaine_=chaine[:i]+str(value)+chaine[i+lon:]
    #print(chaine_)
    return chaine_
    
    return '' 


def get_Par_data(_path):
    D={}
    #params principaux
    D['friction']=''
    D['ngll']=''
    D['fmax']=''
    D['ndof']=''
    D['W']=''
    D['xlim']=''
    D['zlim']=''
    D['nelem']=''
    D['rho']=''
    D['cp']=''
    D['cs']=''
    D['friction']=''
    D['Tn'] =''
    D['Tt']=''
    D['V']=''
    D['kind']=''
    D['Dc']=''
    D['MuS']=''
    D['MuD']=''
    D['a']=''
    D['b']=''
    D['Vstar']=''
    D['Vc']=''
    D['theta']=''
    #param de dimensionnement : ORDER,...
    D['thetaXn']=''
    D['thetaZn']=''
    D['TtXn']=''
    D['TtZn']=''
    D['DcXn']=''
    D['DcZn']=''
    D['aXn']=''
    D['VXn']=''
    entete=["&GENERAL ","&MESH_CART ","&MAT_ELASTIC ","&BC_DYNFLT ","&BC_DYNFLT_RSF ","&BC_DYNFLT_SWF ","&DIST_ORDER0 ","&DIST_PWCONR"]
    def regulariser_elems(var_liste):
        var_liste_=[]
        for e in var_liste:
            ee=e[:]
            if ord(e[0])==9:ee=e[1:]
            if ee[-1]=='\n':var_liste_.append(str(ee[:-1]).replace('d','e'))
            else:var_liste_.append(ee.replace('d','e'))
        return var_liste_
    def remove_elems(splitline):
        splitline_=[]
        for e in splitline:
            if e!='' and e!=' ' and e!='\n':splitline_.append(e)
        return splitline_
    def apply_dist(cnt,f,D,texte,i):
        if cnt>=1:
            var=regulariser_elems(remove_elems(texte[i+2].split(' ')+texte[i+3].split(' ')))
            D[f[0]]= var
            for e in f:
                if e=='theta':
                    D['thetaXn'] = get_value(texte[i+1],'xn')
                if e=='Tt':
                    D['TtXn']    = get_value(texte[i+1],'xn')
                    D['thetaXn'] = get_value(texte[i+1],'xn')
                if e=='V':
                    D['VXn']    = get_value(texte[i+1],'xn')
                    D['thetaXn']= D['VXn']
                    D['TtXn']   = D['VXn'];
                if e=='a':
                    D['aXn']    = get_value(texte[i+1],'xn')
                    D['thetaXn']= D['aXn']
                    D['TtXn']   = D['aXn'];
                if e=='Dc':
                    D['DcXn']    = get_value(texte[i+1],'xn')
                    D['thetaXn'] = D['DcXn']
                    D['TtXn']    = D['DcXn'];
        if cnt>=2:
            var2=texte[i+5].split(' ')+texte[i+6].split(' ')
            var2=remove_elems(var2)
            var2=regulariser_elems(var2)
            D[f[1]]= var2  
        return D
    def count_DIST(ligne,D):#compte le nombre de DIST-ORDER0 SUR LA LIGNE      
        i=0
        cnt=0
        mes_nom_var=[]        
        for ch in ligne:
            if ch=='O' and ligne[i-2]=='=' and ligne[i:i+5]=='ORDER' :
                cnt=cnt+1
                #print(ligne[0:i-1])
                j=3
                nom_var=''
                while ligne[i-j]!=',':
                    nom_var=nom_var+ligne[i-j]
                    j=j+1
                mes_nom_var.append("".join(reversed(nom_var.strip(' ')[1:])) )
            i=i+1
        return cnt,mes_nom_var
    def regulariser(var_):
        var=var_
        if var_[0]=='\\':
            var=var[1:]
        if var_[-1]=='\n' or var_[-1]=='/':
            var=var[:-1]
            if var_[-2]=='\n' or var_[-2]=='/':
                var=var[:-1]
            return var
        return var_       
    texte=[]
    try:
        Par=open(_path+"/"+"savePar.inp","r")
        texte=Par.readlines()
        Par.close()
    except:
        print('!')
        try:
            Par=open(_path+"/"+"Par.inp","r")
            texte=Par.readlines()
            Par.close()
        except:print('!!!erreur lecture Par.inp');return 0         
    i=-1
    for ligne in texte:
        i=i+1#donne acces au numero de ligne
        if ligne[0:9]  ==entete[0]:
            li=ligne[9:].split(',')
            D['ngll']=get_value(ligne,'ngll').replace('d','e')
            D['fmax']=get_value(ligne,'fmax').replace('d','e')
            if len(li)==5:D['ndof']=get_value(ligne,'ndof').replace('d','e')
            elif len(li)==6:
                D['W']   =get_value(ligne,'W').replace('d','e')
                D['ndof']=get_value(ligne,'ndof').replace('d','e')
        elif ligne[0:11] ==entete[1]:
            li=ligne[11:].split(',')
            D['xlim'] = [li[0].strip(' ')[5:].replace('d','e'),li[1].strip(' ').replace('d','e')]
            D['zlim'] = [li[2].strip(' ')[5:].replace('d','e'),li[3].strip(' ').replace('d','e')]
            D['nelem']= [li[4].strip(' ')[6:].replace('d','e'),regulariser(li[5].strip(' ').replace('d','e'))]
        elif ligne[0:13] ==entete[2]:
            li=ligne[13:].split(',')
            D['rho'] = get_value(ligne,'rho').replace('d','e')
            D['cp']  = get_value(ligne,'cp').replace('d','e')
            D['cs']  = get_value(ligne,'cs').replace('d','e')
        elif ligne[0:11]==entete[3]:
            li=ligne[11:].split(',')
            cnt=0
            D['friction']= get_value(ligne,'friction')[1:-1]
            D['Tn']      = get_value(ligne,'Tn').replace('d','e')
            D['Tt']      = regulariser(get_value(ligne,'Tt').replace('d','e'))
            D['V']       = get_value(ligne,'V').replace('d','e')
            cnt,f=count_DIST(ligne,D)
            D=apply_dist(cnt,f,D,texte,i)
        elif ligne[0:15] ==entete[4]:                
            li=ligne[15:].split(',')
            cnt=0#compteur de DIST dans une ligne
            D['kind'] =get_value(ligne,'kind')
            D['Dc']   =get_value(ligne,'Dc').replace('d','e')
            D['MuS']  =get_value(ligne,'MuS').replace('d','e')
            D['a']    =get_value(ligne,'a').replace('d','e')
            D['b']    =get_value(ligne,'b').replace('d','e')
            D['Vstar']=get_value(ligne,'Vstar').replace('d','e')
            D['Vc']   =get_value(ligne,'Vc').replace('d','e')
            D['theta']=regulariser(get_value(ligne,'theta'))
            cnt,f=count_DIST(ligne,D)
            D=apply_dist(cnt,f,D,texte,i)
        elif ligne[0:15] ==entete[5]:
            li=ligne[15:].split(',')
            cnt=0#compteur de DIST dans une ligne
            D['kind'] =get_value(ligne,'kind')
            D['Dc']   =get_value(ligne,'Dc').replace('d','e')
            D['MuS']  =get_value(ligne,'MuS').replace('d','e')
            D['Mud']  =get_value(ligne,'Mud').replace('d','e')
            cnt,f=count_DIST(ligne,D)
            D=apply_dist(cnt,f,D,texte,i)
    #print_(D)
    #print__(D)
    return D

def convertir_Par_data(dic):
    def rm_vide(f):
        while '' in f:
            f.remove('')
        return f
    liste=['ngll','fmax','ndof','W','xlim','zlim','nelem','rho','cp','cs','friction','Vc','Tt','kind','Dc']
    liste=[*liste,*['MuS','MuD','a','b','Vstar','Tn','V','theta','thetaXn','thetaZn','TtXn','TtZn','aXn','VXn'] ]                                       
    for e in liste:
        if dic[e]=='':dic[e]=0
    #Param principaux
    dic['ngll']      = int(dic['ngll'])
    dic['fmax']      = float(dic['fmax'])
    dic['ndof']      = int(dic['ndof'])
    dic['W']         = float(dic['W'])
    try:
        dic['xlim']      = [float(i) for i in dic['xlim']]
        dic['zlim']      = [float(i) for i in dic['zlim']]
        dic['nelem']     = [int(i) for i in dic['nelem']]
    except: print("!mesh error : function convertir_Par_data()")
    dic['rho']       = float(dic['rho'])
    dic['cp']        = float(dic['cp'])
    dic['cs']        = float(dic['cs'])
    dic['friction']  = str(dic['friction'])
    dic['Vc']        = float(dic['Vc'])
    dic['Tt']        = np.float_(rm_vide(dic['Tt']))
    dic['theta']     = np.float_(rm_vide(dic['theta']))
    
    dic['kind']      = int(dic['kind'])
    try:dic['Dc']        = float(dic['Dc'])
    except:
        try:dic['Dc']     = np.float_(rm_vide(dic['Dc']))
        except:ValueError:print('!Dc error : function convertir_Par_data()')
    dic['MuS']       = float(dic['MuS'])
    dic['MuD']       = float(dic['MuD'])
    dic['b']         = float(dic['b'])
    dic['Vstar']     = float(str(dic['Vstar']).replace('d','e'))
    #param de dimensionnement : ORDER,...
    dic['thetaXn']   = int(dic['thetaXn'])
    dic['thetaZn']   = int(dic['thetaZn'])
    dic['TtXn']      = int(dic['TtXn'])
    dic['TtZn']      = int(dic['TtZn'])
    dic['aXn']       = int(dic['aXn'])
    dic['VXn']       = int(dic['VXn'])
    try:dic['Tn']    = float(dic['Tn'])
    except:
       try:dic['Tn'] = np.float_(dic['Tn'])
       except ValueError:print('!Tn error : function convertir_Par_data()')
    try:dic['V']     = float(dic['V'])
    except:
        try:dic['V'] = np.float_(dic['V'])
        except ValueError:print('!V error : function convertir_Par_data()')
    try:dic['a']     =float(dic['a'])
    except:
        try:dic['a'] = np.float_(dic['a'])
        except ValueError:print('!a error : function convertir_Par_data()')
    return dic

def reconvertir_Par_data(dic):
    def ch_to_li(ch):
        if (type(ch)==str):
            tab=((ch.replace('[','')).replace(']',''))
            tab=tab.replace(',',' ')
            tab=tab.replace("'","")
            tab=tab.split(' ')
            if '' in tab:
                tab.remove('')
        else:tab=ch
        return tab
    liste=['ngll','fmax','ndof','W','xlim','zlim','nelem','rho','cp','cs','friction','Vc','Tt','kind','Dc']
    liste=[*liste,*['MuS','MuD','a','b','Vstar','Tn','V','theta','thetaXn','thetaZn','TtXn','TtZn','aXn','VXn'] ] 
    for e in liste:
        dic[e]=((str(dic[e]).replace('e','d')).replace('+','')).strip(' ')

    dic['theta']=ch_to_li(dic['theta'])
    dic['Tt']=ch_to_li(dic['Tt'])
    dic['Dc']=ch_to_li(dic['Dc'])
#    D['ngll']      =str(D['ngll']).replace('e','d')
#    D['fmax']      =str(D['fmax']).replace('e','d')
#    D['ndof']      =str(D['ndof']).replace('e','d')
#    D['W']         =str(D['W']).replace('e','d')
#    D['xlim']      =str(D['xlim']).replace('e','d')
#    D['zlim']      =str(D['zlim']).replace('e','d')
#    D['nelem']     =str(D['nelem']).replace('e','d')
#    D['rho']       =str(D['rho']).replace('e','d')
#    D['cp']        =str(D['cp']).replace('e','d')
#    D['cs']        =str(D['cs']).replace('e','d')
#    D['friction']  =str(D['friction'])
#    D['Tn']        =str(D['Tn']).replace('e','d')
#    #D['Tt']        =[str(f) for f in D['Tt']]
#    D['V']         =str(D['V']).replace('e','d')
#    D['kind']      =str(D['kind']).replace('e','d')
#    D['Dc']        =str(D['Dc']).replace('e','d')
#    D['MuS']       =str(D['MuS']).replace('e','d')
#    D['MuD']       =str(D['MuD']).replace('e','d')
#    D['a']         =str(D['a']).replace('e','d')
#    D['b']         =str(D['b']).replace('e','d')
#    D['Vstar']     =str(D['Vstar']).replace('e','d')
#    D['Vc']        =str(D['Vc']).replace('e','d')
#    #D['theta']     =[str(f) for f in D['theta']]    
#    D['thetaXn']   =str(D['thetaXn'])
#    D['thetaZn']   =str(D['thetaZn'])
#    D['TtXn']      =str(D['TtXn'])
#    D['TtZn']      =str(D['TtZn'])
#    D['aXn']       =str(D['aXn'])
#    D['VXn']       =str(D['VXn'])
    return dic


def export_Par_data(_path,D,ti):
    print('*****************export***************** ')
    global out_path#="D:/STAGEM2/scripts"
    out_dir="newPar"
    entete=["&GENERAL","&MESH_CART","&MAT_ELASTIC","&BC_DYNFLT","&BC_DYNFLT_RSF","&BC_DYNFLT_SWF","&DIST_ORDER0"]
    dic=reconvertir_Par_data(D);print_(dic)
    def li_to_ch(liste):#transf liste en chaine
        chaine=''
        for e in liste:
            chaine=chaine+str(e)+' '
        return chaine   
   
    if out_dir not in [f for f in listdir(out_path) if not isfile(join(out_path, f))]:mkdir(out_path+out_dir)
    N=len([f for f in listdir(out_path+out_dir) if not isfile(join(out_path+out_dir, f))]) 
    tit=str(N)+'_'+ti+str(dic['Dc'])#NOM DU DOSSIER EN SORTIE
    mkdir(out_path+out_dir+'/'+tit)  
    out_=out_path+out_dir+'/'+tit+"/Par.inp"
    texte=[]
    try:#ouverture du template
        Par=open(_path+"/"+"savePar.inp","r")
        texte=Par.readlines()
        Par.close()
    except:
        try:
            Par=open(_path+"/"+"Par.inp","r")
            texte=Par.readlines()
            Par.close()
        except:print('!!!error : lecture fichier');return 0
    out_Par=open(out_,"w")#ouverture du fichier de sortie
    ref=['Tt','Tn','V','theta','Dc','Vstar','MuS','a','b']
    i=0
    while i<len(texte):
        txt=texte[i]
        dic['TtXn']=int(dic['TtXn'])
        dic['thetaXn']=int(dic['thetaXn'])
        if 'Tt' in ref and texte[i].split(' ')[0]==entete[3]:#&BC_DYNFLT 
            if 'Tn' in ref:txt=update_texte(texte[i],'Tn',dic['Tn'])
            if 'V' in ref: txt=update_texte(txt,'V',dic['V'])
            out_Par.write(txt)
            txt_1=update_texte(texte[i+1],'xn',dic['TtXn'])
            out_Par.write(txt_1)
            #D['TtZn']=int(get_value(texte[i+1],'xn'))
            #zn=int(get_value(texte[i+1],'TtZn'))
            out_Par.write(li_to_ch(dic["Tt"][0:dic['TtXn']-1]))
            out_Par.write('\n')
            out_Par.write(li_to_ch(dic["Tt"][dic['TtXn']-1:]))
            out_Par.write('\n')
            i=i+3
        elif 'theta' in ref and texte[i].split(' ')[0]==entete[4]:#&BC_DYNFLT_RSF
            if 'Vstar' in ref:txt=update_texte(txt,'Vstar',dic['Vstar'])
            if 'MuS' in ref:txt=update_texte(txt,'MuS',dic['MuS'])
            if 'a' in ref:txt=update_texte(txt,'a',dic['a'])
            if 'b' in ref:txt=update_texte(txt,'b',dic['b']) 
            out_Par.write(txt)
            txt_1=update_texte(texte[i+1],'xn',dic['thetaXn'])
            if type(D['Dc'])!=list:
                if 'Dc' in ref:txt=update_texte(texte[i],'Dc',dic['Dc'])                  
                out_Par.write(txt_1)
                out_Par.write(li_to_ch(dic["theta"][0:dic['thetaXn']-1]))
                out_Par.write('\n')
                out_Par.write(li_to_ch(dic["theta"][dic['thetaXn']-1:]))
                out_Par.write('\n')
                i=i+3
            else:# Dc est une distribution
                if len(D['Dc'])!=1:#si Dc n  est pas un scalaire (securite )
                    out_Par.write(txt_1)
                    out_Par.write(li_to_ch(dic["Dc"][0:dic['thetaXn']-1]))
                    out_Par.write('\n')
                    out_Par.write(li_to_ch(dic["Dc"][dic['thetaXn']-1:]))
                    out_Par.write('\n')
                    i=i+3
                out_Par.write(txt_1)
                out_Par.write(li_to_ch(dic["theta"][0:dic['thetaXn']-1]))
                out_Par.write('\n')
                out_Par.write(li_to_ch(dic["theta"][dic['thetaXn']-1:]))
                out_Par.write('\n')
                i=i+3
                
        else:out_Par.write(texte[i])
        i=i+1
    out_Par.close()
    return out_,tit


#-----------------------------------------------------------------------------------
#                                    '''CALCULS'''
#-----------------------------------------------------------------------------------  


def steady_state(dic):
    s_s=[]
    if len(dic["theta"])==1:s_s=dic["theta"][0]*dic["Vstar"]/dic["Dc"]
    else:#si teta vecteur
        if type(dic['Dc'])!=list and  type(dic['Dc'])!=np.ndarray :#si Dc est scalaire
            for e in dic["theta"][dic['thetaXn']-1:]:
                #print(type(e),type(D["Dc"]),type(D["Vstar"]))
                s_s.append(float(e)*float(dic["V"])/float(dic["Dc"]))
        elif len(dic["theta"])==len(dic["Dc"]):#theta vecteur , Dc vecteur , de meme taille
            for e,d in zip(dic["theta"][dic['thetaXn']-1:],dic["Dc"][dic['thetaXn']-1:]):
                #print(type(e),type(D["Dc"]),type(D["Vstar"]))
                s_s.append(float(e)*float(dic["V"])/float(d))
    return s_s


def rsf_mu(dic):
    #Tt=Tn*mu
    #
    mu    =0
    kind  =(dic["kind"])
    v     =(dic["V"])
    Dc    =(dic["Dc"])
    theta =(dic["theta"])
    mus   =(dic["MuS"])
    a     =(dic["a"])
    b     =(dic["b"])
    Vstar =(dic["Vstar"])
    Vc    =(dic["Vc"])
    tn    =float((dic["Tn"]))
    if type(dic["theta"])!=list:#si theta is scalaire
        theta=dic["theta"][0]
        if kind==1:mu= mus +a*abs(v)/(abs(v)+Vstar) - b*theta/(theta+Dc)#Strong velocity-weakening at high speed  as in Ampuero and Ben-Zion        
        elif kind==2 or kind==3:mu= a*np.arcsinh(abs(v)/(2*Vstar)*np.exp((mus+b*np.log(Vstar*theta/Dc))/a))#Kaneko et al. (2008) Eq. 15 (regularize at v=0 as per Laupsta et al. (2000))
        elif kind==4:mu= a*np.asinh(abs(v)/(2*Vstar)*np.exp((mus+b*np.log(Vc*theta/Dc+1))/a))
        dic['Tt']= -tn*mu
    else:#theta is vecteur
        if type(Dc)!=list:#si Dc scalaire,theta is vecteur
            i=-1
            theta=dic["theta"][dic["thetaXn"]-1:]
            mu=[0]*(len(theta))
            tt0=[0]*(len(theta))
            dic['Tt']=[0]*(len(theta)+dic['thetaXn']-1)
            dic['Tt'][0:dic['thetaXn']]=dic['theta'][0:dic['thetaXn']]
            for t in theta:
                i=i+1
                if kind==1:mu[i]= mus +a*abs(v)/(abs(v)+Vstar) - b*t/(t+Dc)
                elif kind==2 or kind==3:
                    mu[i]= a*np.arcsinh(abs(v)/(2*Vstar)*np.exp((mus+b*np.log(Vstar*t/Dc))/a))
                elif kind==4:mu[i]= a*np.arcsinh(abs(v)/(2*Vstar)*np.exp((mus+b*np.log(Vc*t/Dc+1))/a))
                tt0[i]=-tn*mu[i]
            dic['Tt'][dic["TtXn"]-1:]= tt0
        elif len(Dc)==len(theta):#si Dc ,theta is vecteur et sont de meme taille
            i=-1
            theta=dic["theta"][dic["thetaXn"]-1:]
            Dc=dic["Dc"][dic["thetaXn"]-1:]#on definit en plus les Dci
            mu=[0]*(len(theta))
            tt0=[0]*(len(theta))
            dic['Tt']=[0]*(len(theta)+dic['thetaXn']-1)
            dic['Tt'][0:dic['thetaXn']]=dic['theta'][0:dic['thetaXn']]
            for t,dc in zip(theta,Dc):
                i=i+1
                if kind==1:mu[i]= mus +a*abs(v)/(abs(v)+Vstar) - b*t/(t+dc)
                elif kind==2 or kind==3:mu[i]= a*np.arcsinh(abs(v)/(2*Vstar)*np.exp((mus+b*np.log(Vstar*t/dc))/a))
                elif kind==4:mu[i]= a*np.arcsinh(abs(v)/(2*Vstar)*np.exp((mus+b*np.log(Vc*t/dc+1))/a))
                tt0[i]=-tn*mu[i]
            dic['Tt'][dic["TtXn"]-1:]= tt0
        else:print('problemme rsf_mu')
    return dic

def swf_mu(f,kind):
    mu=0
    mus=f["mus"]
    mud=f["mud"]
    alpha=f["alpha"]
    theta=f["theta"]
    dc=f["Dc"]
    if (kind==1): #!-- linear slip weakening:
        mu = mus -(mus-mud)*np.min(theta/dc)
    else: #!-- exponential slip weakening:
        mu = mud -(mud-mus)*np.exp(-theta/dc)
    mu = mu + alpha * theta
    return mu








