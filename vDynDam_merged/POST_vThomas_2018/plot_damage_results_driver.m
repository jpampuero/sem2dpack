% ========================================================================
% Script to plot the outputs from a dynamic damage material, SEM2DPAK, v18

% Marion Thomas, last modified August 2018

%!! you need to have the GHOStsCRIPT (http://pages.uoregon.edu/koch)
%installed on your computer and pdftops packages installed on your computer
% ========================================================================

% It compares the results from the damage simulations with those of an
% elastic case for :

% 1. Damage with profiles across the fault and log plot
% 2. Snapshots (D,cs,cp,etc..) with cumulative slip and slip velocity (FIGURES PAPER MONO)
% 3. cs/cp reduction in % (FIGURES PAPER GJI)
% 4. Particules velocity (snpashot)
% 5. strain, stress (snpashot)
% 6. Invariants (snpashot)
% 7. Multiple snapshots
    % a. Damage
    % b. cs*
    % c. cp*
    % d. dldt
    % e. Regime
    % f. velocity
% 8. SEISMOMETERS - velocity
    % a. particules velocity
    % b. Fourrier analysis (FAS)
    % c. FAS ratio
% 9. SEISMOMETERS - velocity
    % a. particules velocity
    % b. Fourrier analysis (FAS)
    % c. FAS ratio
% 9. BEHAVIOR ON THE FAULT (slip, slip rate, normal and shear stress)

%Variables that can be plotted as snapshot
    % D:           damage density
    % cs:          S-waves velocity
    % cpp:         P-waves velocity
    % R:           Regime
    % ax,az:       particules acceleration
    % vx,vz:       particules velocity
    % nv:          norm of  particules velocity
    % e11,e22,e12: strain field
    % s11,s22,s12: stress field
    % invI:        first invariant of the stress field
    % invII:       second invariant of the stress field
    % CFF:         Coulomb stress field (invII+0.6*invI)/Drucker-Prager criterion

% CALLS: getsimulinfo.m; gettimesteps.m; sem2d_read_specgrid.m; sem2d_read_fault.m sem2d_read_seis.m;
%sem2d_spapshot_read.m,plot_snapshots_dyD;plot_snapshots_cscp;
%sem2d_snapshot_plot_km; plot_snapshots; plot_snapshots_dyD_multistep; plot_snapshots_sliprate_multistep
%plot_snapshots_dyD_movie.m


%==========================================================================
clc;clear;close all; 

%% ADD DIRECTORIES TO PATH

rootD = [pwd '/'];
add_subdir_to_path(rootD)

%%
% ========================== to be modified ==============================
%                        PARAMETERS FOR THE MODEL
% ========================================================================

%% PATHS
%directories were the model ouptut files are
dir0= '/Users/florescuba/Documents/SEM2D/vDynDam_merged/EXAMPLES/Dyn_Damage/res30m/';

%Sub-folder
% nfault='roughF1e_2/';
% nfault='roughF5e_3/';
% nfault='roughF1e_3/';
nfault='flatF_ang0/';
dir1=[dir0,nfault];

%Roughness
%alpha = 1e-2 : 2043060600 ; 2043063178 ; 2043093826 ; 2043104396 ; 2043231839
%alpha = 5e-3 : 2043194986 ; 2043197257 ; 2043199051 ; 2043200602 ; 2043206198
%alpha = 1e-3 : 2043156479 ; 2043159460 ; 2043161324 ; 2043169931 ; 2043177427
% seed = 113160203;

%Damage
% datadirD = [dir1 'dydmg_roughF_' num2str(seed) '_1D0/'];
datadirD = [dir1 'dydmg_2D0/'];

%Elasticity
% datadirE = [dir1 'nodmg_roughF_' num2str(seed) '/'];
datadirE = [dir1 'nodmg/'];

%Flat fault elasticity
% n3='flatF_ang0/';
% datadirF = [dir0,n3,'nodmg/'];
datadirF =datadirE;

%% VARIABLES TO BE DEFINED

xb=[-6e3,6e3];      %x boundaries for the plot
zb=[-1.0e3,1.0e3];  %z boundaries for the plot
save_tag = 0;       %save figure (1=save, 0=just plot, 2=save just .fig)
tplot = 2.7;        %Max time (in s) at which we want to plot the variable   
tplotF = 2.7;       %Max time (in s) at which we want to plot the variable for the flat case
tstart = 2.0;       %starting time to plot the variables along the fault
ndt = 5;            %number of iscochrone for the slip and slip velocity plots
fs = 14;            %font size;

%colormaps
load('stresscmap','mycmapS');NmycmapS=ones(1,3)*size(mycmapS,1);
load('parula_extended','mycmap');Nmycmap=ones(1,3)*size(mycmap,1);
load('vik_25');vik_25(13,:) = repmat([1 1 1],1,1);Nvik_25=ones(1,3)*size(vik_25,1);
cmapD = CubeHelix(256,0.4,-1.5,2,0.7);NcmapD=ones(1,3)*size(cmapD,1);
cmapreg=[1 1 1;0 0 0; 1 0.9 0];Ncmapreg=ones(1,3)*size(cmapreg,1);
cmapreg2=[1 1 1; 0 0 0; 0.5 0.5 0.5];Ncmapreg2=ones(1,3)*size(cmapreg2,1);
mycmap3=[flip(mycmap(:,1)),flip(mycmap(:,2)),flip(mycmap(:,2))];Nmycmap3=ones(1,3)*size(mycmap3,1);
cmapcscp=[1,1,1;flip(mycmap)];Nmycmap4=ones(1,3)*size(cmapcscp,1);
cgray = flip(gray);Ncgray = ones(1,3)*size(cgray,1);

%% 
% =========================================================================
%                           PLOTTING FUNCTIONS
% =========================================================================
%
%% DATA DOWNLOAD 

%Info for plots 
load_info

%Load data from 
load_data

%ouput names
output_name

%Check if the user want to continue
if checkvarD==1
    disp(' '); disp(['it might replace existing files in ', namefoldD])
	m=input('Do you want to continue? y/n:','s');
	if strcmp(m,'N')==1 || strcmp(m,'n')==1; return; end
end

disp(' ')
disp('========================================');
disp(['Figures for : ',strtok(datadirD(length(dir1)+1:end-1))]);
disp('========================================');

pro_zone=1e3;%1.2589e+03

%% 1. Snaphshot of Damage with profiles across the fault and log plot

%options:
norma=pro_zone;             %normalization for the x- and z-axis
lim_axis_var = [xb,zb,1]; %min/max for x-axis, z-axis and the aspect ratio between axis (x/z)
limD = [0.1 1.0];           %minimum and maximum value to be plotted
np=20;                      %number of profiles
smooth_op=0;                %smoothing option for snapshot: 1 for yes, 0 for no.

%profiles
options=[lim_axis_var,limD,fs,infoP,np,norma,tplot,smooth_op];
plot_profil_dyD(datadirD,gridD,ooD,TD,faultD,save_tag,namefoldD,options,'white',cmapD)

%% 2. Snapshots (D,cs,cp,etc..) with cumulative slip and slip velocity (FIGURES PAPER MONO)

% %options:
% norma=pro_zone;             %normalization for the axis
% lim_axis_var = [xb,zb,0.5]; %min/max for x-axis, y-axis and the aspect ratio between axis (x/y)
% lim_var = [0.1 1.0];        %minimum and maximum value to be plotted
% tstp=0.35;                  %time steps for the slip and slip velocity plots
% ym_d = 10;                  %max slip for slip  plot
% ym_v = 12;                  %max slip for slip rate plot
% fs = 14;                    %font size;
% as2 = 2.16/(((xb(2)-xb(1))/1e3)/ym_d); %1.2;%aspect ratio between axis (x/y) for slip
% as3 = 2.16/(((xb(2)-xb(1))/1e3)/ym_v); %1.8;%aspect ratio between axis (x/y) for slip rate
% 
% %combine options
% options=[lim_axis_var,lim_var,tstp,ym_d,ym_v,fs,infoP,as2,as3,norma,tplot];
% 
% % plot: you can plot D,cs,cp,dldt,KI,A,B,C,R
% plot_snapshots_dyD(datadirD,gridD,ooD,'D',TD,faultD,faultE,save_tag,namefoldD,options,'white',cmapD)

%% 3. cs/cp reduction in % (FIGURES PAPER GJI)

%options:
norma=pro_zone;             %normalization for the axis
lim_axis_var = [xb,zb,0.5]; %min/max for x-axis, y-axis and the aspect ratio between axis (x/y)
lim_var1 = [0 35];          %minimum and maximum value to be plotted
lim_var2 = [0 35];          %minimum and maximum value to be plotted
DThres = 0.3;               %threshold to plot cs/cp variation
smooth_op=1;                %smoothing option for snapshot: 0 for NO, 1 for slight smoothing,
                            %2 for high smoothing like Kurama

%combine options
options=[lim_axis_var,fs,infoP,norma,tplot,DThres,smooth_op,lim_var1,lim_var2];
% options=[lim_axis_var,fs,infoP,norma,tplot,DThres,smooth_op];

%plot
plot_snapshots_cscp(datadirD,gridD,ooD,TD,faultD,save_tag,namefoldD,options,'black',cmapcscp)

%% 4. Particules velocity

%options:
norma=pro_zone;             %normalization for the axis
lim_axis_var = [xb,zb,0.5];   %min/max for x-axis, y-axis and the aspect ratio between axis (x/y)
lim_var = [-3 3];           %minimum and maximum value to be plotted
ym_v = 15;                  %max value for slip rate plot 
smooth_op=0;                %smoothing option for snapshot: 0 for NO, 1 for slight smoothing,
                            %2 for high smoothing like Kurama

%plot
options=[lim_axis_var,fs,infoP,ym_v,norma,2.3,smooth_op,lim_var,lim_var,[0 3]];%,lim_var
clear var;
var{1}='vx';var{2}='vz';var{3}='nv';
% cbar=[mycmap;mycmap;mycmap;Nmycmap;Nmycmap;Nmycmap];
cbar=[vik_25;vik_25;vik_25;Nvik_25;Nvik_25;Nvik_25];
plot_snapshots_nvar(datadirD,gridD,ooD,var,TD,faultD,save_tag,namefoldD,options,cbar,'black')
plot_snapshots_nvar(datadirE,gridE,ooE,var,TE,faultE,save_tag,namefoldE,options,cbar,'black')

%% 5. strain, stress

% %options:
% norma=pro_zone;             %normalization for the axis
% lim_axis_var = [xb,zb,0.5]; %min/max for x-axis, y-axis and the aspect ratio between axis (x/y)
% lim_var = [-1.2e-3 1.2e-3];       %minimum and maximum value to be plotted (in MPa for stress and strain)
% ym_v = 15;                  %max value for slip rate plot  
% smooth_op=0;                %smoothing option for snapshot: 0 for NO, 1 for slight smoothing,
%                             %2 for high smoothing like Kurama
% 
% % %a. strain
% options=[lim_axis_var,fs,infoP,ym_v,norma,tplot,smooth_op,lim_var,lim_var,[-5e-3 5e-3]];
% clear var;
% var{1}='e22';var{2}='e12';var{3}='e11';
% cbar=[mycmap;mycmap;mycmap;Nmycmap;Nmycmap;Nmycmap];
% plot_snapshots_nvar(datadirD,gridD,ooD,var,TD,faultD,save_tag,namefoldD,options,cbar,'white')
% plot_snapshots_nvarE(datadirE,gridE,ooE,var,TE,faultE,save_tag,namefoldE,options,cbar,'white')
% 
% %b stress
% options=[lim_axis_var,fs,infoP,ym_v,norma,tplot,smooth_op,lim_var,lim_var,lim_var];
% clear var;
% var{1}='s22';var{2}='s12';var{3}='s11';
% cbar=[mycmap;mycmap;mycmap;Nmycmap;Nmycmap;Nmycmap];
% plot_snapshots_nvar(datadirD,gridD,ooD,var,TD,faultD,save_tag,namefoldD,options,cbar,'white')
% plot_snapshots_nvarE(datadirE,gridE,ooE,var,TE,faultE,save_tag,namefoldE,options,cbar,'white')

%% 6. Invariants, drucker-praguer and damage

%options:
norma=pro_zone;             %normalization for the axis
lim_axis_var = [xb,zb,0.5]; %min/max for x-axis, z-axis and the aspect ratio between axis (x/z)
ym_v = 15;                  %max value for slip rate plot 
t1=tplot;                     %Time at which we want to plot the variable
smooth_op=2;                %smoothing option for snapshot: 0 for NO, 1 for slight smoothing,
                            %2 for high smoothing like Kurama

%Colors bars values; 
%If you don't set any in the "options" variable, it will take automatically the min/max values for the variable to be plotted
lim_varI = [-1 1];          %min/max values to be ploted for 1st invariant
lim_varII = [-1.7 1.7];     %min/max values to be ploted for 2nd invariant
lim_CFF = [-50 50];     	%min/max values to be ploted for drucker-praguer criterion
lim_R = [0.5 3.5];          %min/max values to be ploted for regime
lim_D = [0.2 1];            %min/max values to be ploted for cumulative damage at the time of the snapshot

%plot drucker-praguer criterion
options=[lim_axis_var,fs,infoP,ym_v,norma,t1,smooth_op,lim_CFF];
clear var;
var{1}='CFF';% var{1}='D';var{2}='R';
cbar=[vik_25;Nvik_25];
plot_snapshots_nvar(datadirD,gridD,ooD,var,TD,faultD,save_tag,namefoldD,options,cbar,'black')

% %plot cumulative damage and damage at the time of the snapshot (ie. Regime)
% smooth_op=0;
% options=[lim_axis_var,fs,infoP,ym_v,norma,t1,smooth_op,lim_D,lim_R];
% clear var;
% var{1}='D';var{2}='R';
% cbar=[cgray;cmapreg;Ncgray;Ncmapreg];
% plot_snapshots_nvar(datadirD,gridD,ooD,var,TD,faultD,save_tag,namefoldD,options,cbar,'black')

% % plot first and second invariants
% smooth_op=0;
% options=[lim_axis_var,fs,infoP,ym_v,norma,t1,smooth_op,lim_varI,lim_varII];%,lim_CFF];
% clear var; 
% var{1}='invI';var{2}='invII';%var{3}='CFF';
% cbar=[vik_25;vik_25;vik_25;Nvik_25;Nvik_25;Nvik_25];
% plot_snapshots_nvar(datadirD,gridD,ooD,var,TD,faultD,save_tag,namefoldD,options,cbar,'black')
% plot_snapshots_nvarE(datadirE,gridE,ooE,var,TE,faultE,save_tag,namefoldE,options,cbar,'black')
% 
%% 7. Multiple snapshots

%options:
n=3;                        %number of subplot
norma=pro_zone;             %normalization for the axis
lim_axis_var = [xb,[-2e3,2e3],1]; %min/max for x-axis, y-axis and the aspect ratio between axis (x/y)
ym_v = 15;                  %max slip for slip rate plot
t1=tstart; t2=tplot;             %starting and ending time for snapshots
smooth_op=0;                %smoothing option for snapshot: 0 for NO, 1 for slight smoothing,
                            %2 for high smoothing like Kurama

% a. Damage
lim_D = [0.1 1.0];        %D: minimum and maximum value to be plotted
options=[lim_axis_var,fs,ym_v,norma,t1,t2,smooth_op,lim_D];
plot_snapshots_dyD_multistep(datadirD,gridD,n,ooD,'D',TD,faultD,save_tag,namefoldD,options,cmapD,'white')

% % b. cs*
% lim_var = [0 40];         %minimum and maximum value to be plotted
% options=[lim_axis_var,fs,ym_v,norma,tstart,tplot,smooth_op,lim_var];
% plot_snapshots_dyD_multistep(datadirD,gridD,n,ooD,'cs',TD,faultD,save_tag,namefoldD,options,cmapcscp,'black')
% 
% % c. cp*
% lim_var = [0 40];         %minimum and maximum value to be plotted
% options=[lim_axis_var,fs,ym_v,norma,tstart,tplot,smooth_op,lim_var];
% plot_snapshots_dyD_multistep(datadirD,gridD,n,ooD,'cp',TD,faultD,save_tag,namefoldD,options,cmapcscp,'black')
% 
% % d. Regime
% lim_var = [0.5 3.5];      %minimum and maximum value to be ploted
% options=[lim_axis_var,fs,ym_v,norma,t1,t2,smooth_op,lim_var];
% plot_snapshots_dyD_multistep(datadirD,gridD,n,ooD,'R',TD,faultD,save_tag,namefoldD,options,cmapreg,'white')
% 
% % e. particules velocity
% lim_var = [0 9];          %minimum and maximum value to be ploted
% options=[lim_axis_var,fs,ym_v,norma,t1,t2,smooth_op,lim_var];
% plot_snapshots_dyD_multistep(datadirD,gridD,n,ooD,'nv',TD,faultD,save_tag,namefoldD,options,mycmap,'white')
% 
% % f. CFF
% lim_var = [-40 40];       %minimum and maximum value to be ploted
% options=[lim_axis_var,fs,ym_v,norma,t1,t2,lim_var];
% plot_snapshots_dyD_multistep(datadirD,gridD,n,ooD,'CFF',TD,faultD,save_tag,namefoldD,options,vik_25,'black')

%% 8. SEISMOMETERS - velocity and Fourrier analysis
% clc; close all

% Points history coordinates, (x,z) (should corresponds to the input file sensors.ini)
Pcoord=[-7500 125;-7500 500;-7500 4000;7500 -125; 7500 -500;7500 -40500;...
    -5000 125;-5000 500;-5000 4000; 5000 -125; 5000 -500;5000 -4000]; %FF
% Pcoord=[5000 -125; 5000 -500;5000 -2500];%-5000 125;-5000 500;-5000 2500]; %FF
% Pcoord=[7500 -125; 7500 -500;7500 -2500;5000 -125; 5000 -500;5000 -2500];%-5000 125;-5000 500;-5000 2500]; %FF



% %data
% data=[sismo;sismoE;sismoF];             %sismo data for the simulations you want to plot
% fault=[faultD;faultE;faultF]; 
% tend=[simuD;simuE;simuF];
% % tend=ones(3,1)*min([simuD;simuE;simuF]);
data=[sismo;sismoE];             %sismo data for the simulations you want to plot
fault=[faultD;faultE]; 
tend=[simuD;simuE];

%options I:
norma=pro_zone;             %normalization for the axis
fs = 14;                    %font size;

% a. plot velocity
field = 'z';                %fields to plot (x,z or magnitude 'm')
lim_varX = [1.0 3.5];     %min/max time for the seismograms
lim_varY = [-5 5];          %min/max amplitude for the seismograms 
options=[fs,norma,3,lim_varY,lim_varX];
plot_velocity(Pcoord,data,fault,tend,save_tag,namefoldD,options,field)

% b. plot FAS
field = 'z';   
fcutoff = 350;              %Frequency cut off
taperlengthpercent = 5;     %Cosine Taper`
lim_varY = [-5 1];          %min/max magnitude of the fourier transform
options=[fs,norma,fcutoff,taperlengthpercent,lim_varY];
plot_fourrier_analysis(Pcoord,data,fault,tend,save_tag,namefoldD,options,field)

% c. plot FAS ratio
fcutoff = 350;              %Frequency cut off
taperlengthpercent = 5;     %Cosine Taper`
lim_varY = [-3 3];          %min/max magnitude of the fourier transform
options=[fs,norma,fcutoff,taperlengthpercent,lim_varY];
plot_FAS_ratio(Pcoord,data,fault,tend,save_tag,namefoldD,options,field)

% % d. plot fault-normal and fault-parallel velocity
% lim_varX = [0.5 simuD];	%min/max time for the seismograms
% lim_varY = [-7 7];      %min/max amplitude for the seismograms 
% options=[fs,norma,simuD];%,lim_varY,lim_varX];
% plot_velocity_VpVn(Pcoord,[sismo;sismoE],[faultD;faultE],[simuD;simuE],save_tag,namefoldD,options)% 
% 
% e. plot FAS for fault-normal and fault-parallel velocity
field = 'z';   
fcutoff = 350;              %Frequency cut off
taperlengthpercent = 5;     %Cosine Taper`
lim_varY = [-5 1];          %min/max magnitude of the fourier transform
options=[fs,norma,fcutoff,taperlengthpercent,lim_varY];
% plot_fourrier_analysis_VpVn(Pcoord,[sismo;sismoE],[faultD;faultE],[simuT;simuT],save_tag,namefoldD,options,field)
plot_fourrier_analysis_VpVn(Pcoord,data,fault,tend,save_tag,namefoldD,options,field)
% 
% % c. plot FAS ratio
% fcutoff = 300;              %Frequency cut off
% taperlengthpercent = 5;     %Cosine Taper`
% lim_varY = [-3 3];          %min/max magnitude of the fourier transform
% options=[fs,norma,fcutoff,taperlengthpercent,lim_varY];
% plot_FAS_ratio_VpVn(Pcoord,data,fault,tend,save_tag,namefoldD,options)
% 
% sismo_location
norma=pro_zone;             %normalization for the axis
lim_axis_var = [ [-8e3 8e3],[-3960,4080],1];   %min/max for x-axis, y-axis and the aspect ratio between axis (x/y)
options=[lim_axis_var,fs,infoP,norma,lim_varX];
sismo_location(datadirD,gridD,ooD,TD,Pcoord,data,fault,tend,save_tag,namefoldD,options,'white',cmapD)
 
%% 9. BEHAVIOR ON THE FAULT

%options:
ndt = 8;                %number of iscochrone
norma=pro_zone;         %normalization for the axis
ym_d = 10;              %max slip for slip  plot 
ym_v = 15;              %max slip for slip rate plot 
ym_sn = [-1.3 1.3];     %min/max normalized normal stress
ym_st = [-3.6 3.0];     %min/max normalized shear stress

%combine options
options=[xb,ndt,ym_d,ym_v,ym_sn,ym_st,fs,norma,0,simuD];
plot_fault_parameters(faultD,faultE,save_tag,namefoldD,options)

%% 10. MOVIES (D and slip rate on the fault)

% %Create folders.
% mkdir([root_directory,n1,name1a,'/movie/'])
% mkdir([root_directory,n1,name1a,'/movie/D/'])
% 
% % Damage
% %options:
% xb=[-9e3,9e3];
% yb=[-1e3,1e3];
% lim_axis_var = [xb,yb,0.25]; %min/max for x-axis, y-axis and the aspect ratio between axis (x/y)
% as2 = 1.2;                  %aspect ratio between axis (x/y) for snapshot
% lim_var = [0.0 1];          %minimum and maximum value to be ploted
% ym_v = 14;                  %max slip for slip rate plot
% fs = 10;                    %font size;
% co1='Dynamic Damage'; clr1='red';
% co2='Linear elastic'; clr2='red';
% %combine options
% options=[lim_axis_var,as2,lim_var,ym_v,fs,info1,info2,tplot];
% % MOVIE snapshots
% plot_snapshots_dyD_movie(datadir1,grid1,oo1,'D',T1,name1,data1,time1,namef,options)

% %Velocity field
% mkdir([root_directory,n1,name1a,'/movie/v/'])
% co1='Dynamic Damage evolution'; clr1='red';
% co2='No Damage evolution'; clr2='white';
% plot_snapshots_movie(datadir1,grid1,numel(oo1),oo1,'v',3,0,T1,name1,xb,yb,data1,time1,co1,co2,clr1,clr2,n1)

%% POINTS OUTPUTS DOWNLOADS from snapshots

% %Download seismograms from snapshots
% aa = read_points(datadir1,grid1,oo1,'a',T1,Pcoord);
% aaE = read_points(datadir2,grid2,oo2,'a',T2,Pcoord);
% vv = read_points(datadir1,grid1,oo1,'v',T1,Pcoord);
% vvE = read_points(datadir2,grid2,oo2,'v',T2,Pcoord);
% 
% %options:
% xb=[-9e3,9e3];
% yb=[-2e3,2e3];
% lim_axis_var = [xb,yb,0.25]; %min/max for x-axis, y-axis and the aspect ratio between axis (x/y)
% lim_var = [0.0 1];         %minimum and maximum value to be ploted
% fs = 14;                    %font size;
% 
% %options:
% field = 'z';    % Fields to plot (x,z or magnitude 'm')
% options=[lim_axis_var,lim_var,fs,info1,info2];
% %plots
% plot_velocity_from_snapshot(datadir1,grid1,oo1,Pcoord,sismo,sismoE,timesis,timesisE,name1,save_tag,namef,options,field)
% 
% %options:
% field = 'x';    % Fields to plot (x,z or magnitude 'm')
% options=[lim_axis_var,lim_var,fs,info1,info2];
% %plots
% plot_velocity_from_snapshot(datadir1,grid1,oo1,Pcoord,sismo,sismoE,timesis,timesisE,name1,save_tag,namef,options,field)

