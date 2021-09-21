% ========================================================================
% Script to build the grid and create a SEM2DPACK input file for a rough or 
% a straight fault, along the x-axis, with dynamic damage. It can create 
% longitudinal variation of bulk properties.

% Marion Thomas, last modified August 2018

% CALLS: add_subdir_to_path; compute_ini_param; check_D0 ; compute_ini_CS_CP
%        fault_gen; mesh2d_generate; output_name; mesh2d_write; dtau_comp
%        create_inputfile
%==========================================================================
clc;clear;close all; 

%% ADD DIRECTORIES TO PATH

rootD = [pwd '/'];
add_subdir_to_path(rootD)

%%
% ========================== to be modified ==============================
%                        PARAMETERS FOR THE MODEL
% ========================================================================

%% MESH PARAMETERS

%Mesh grid
BOX = [-32e3 32e3 -3.0e3 3.0e3];  % domain limits in m [xmin xmax zmin zmax]
h = 30;                         % element size in m (resolution)
Q = 4;                          % quad element type, identified by its number of nodes: 4 or 9.
BOX=resocorrec(BOX,h,Q);        % correct the dimension of the domain based on the resolution

%simulations parameters
ngll=3;     % Number of GLL nodes per edge on each spectral element
fmax=4d0;   % maximum frequency to be resolved
ndof=2;     % Number of degrees of freedom per node: 1 for antiplane, 2 for inplane

%% FAULT GEOMETRY

%Flat fault
% Fl = [BOX(1)+3e3 0]; %(x,z) Coordinates of the fault on the left side
% Fr = [BOX(2)-2e3 0]; %(x,z) Coordinates of the fault on the right side
Fl = [BOX(1) 0]; %(x,z) Coordinates of the fault on the left side
Fr = [BOX(2) 0]; %(x,z) Coordinates of the fault on the right side

% rought fault
% rms=5e-3;                   % root mean square (rms) height fluctuations of order "rms" times the profile length.
% nfault='roughF5e_3';    	% name to assess the roughness (for the mesh grid name)
% seed=[23059117];                    % Give a number or leave it empty for random generation of fault roughness
% Seed used
% alpha = 1e-2 : 2043231839 ; 2043060600 ; 2043063178 ; 2043093826 ; 2043104396
% alpha = 5e-3 : 2043194986 ; 2043197257 ; 2043199051 ; 2043200602 ; 2043206198
% alpha = 1e-3 : 2043177427 ; 2043156479 ; 2043159460 ; 2043161324 ; 2043169931 

%% FAULT PARAMETERS

%Slip weakening law
MuS=0.6d0;          % static coefficient of friction
MuD=0.1d0;          % dynamic coefficient of friction
Dc=1d0;

%Normal stress response for dynamic faults (cf. SEM2DPACK users'guide)
kind=2;     % 2 = Prakash-Clifton with regularizing time scale
T=4.0e-3;   % Regularization time scale if kind=2

%nucleation prone parch
nucsize = 1.5e3;            %radius (in m) of the prone patch (around "nucloca")
% nuclocX = -25e3;              %give the location of the prone patch along the x axis
nuclocX = 0e3;              %give the location of the prone patch along the x axis
correctionfactor = 1.001;   %correction factor to MuS to compute the assign the
                            %shear stress so that it is above the fault strengt

% To compute initial strain, SXZ, SXX AND MU0
PS=[-8;-5;-3.4]*1e7; % Principal stress
ang=30;              % based on Mohr coulomb circle. Angle of the normal to the fault with 
                     % principles stress (sigma 1) has to be greater than 45

                     
%% BULK PROPERTIES & DAMAGE LAW

%Bulk properties: you may assign 2, if you want to have a different
%material for the bottom (first value) and the top (2nd value)
% rho=[2700,3000];           % density
% cs=[3.11508e3,3.25e3];       % s-waves velocity
% cp=[5.6e3,5.84e3];           % p-waves velocity
rho=[2700];           % density
cs=[3.12e3];       % s-waves velocity
cp=[5.6e3];           % p-waves velocity

%Properties for Damage law
a0=60;              % initial radius of the crack (in m)
fs=MuS;             % coefficient of friction
vm = 1.5832e3;      % Branching speed
beta=0.1e0; 
omega=2.0e0;
KicSS=1.2d6;        % quasi-static fracture toughness (MPa x m1/2)
fitKicD = 5e9;

%% INITIAL DAMAGE

%Default/background values for initial Damage (D0) on both side of the fault
D0top=0.1;
D0bot=0.1;

% if you want to vary D0 inside the medium, there are few possibilities. Set
% the tag to 1 if you want to use one options. With this code, the variations are 
% always perpendicular to the fault.
D0p=[];D0pz=[]; %Default values


% YOUR OWN PECULIAR distribution of D0
%-------------------------------------
% 1)Define a D0 profile across the fault (ownD0). 
% 2) Define ownD0p: it corresponds to the max distance away from the fault for 
%    which ownD0(i) is applicable : ownD0p<0 = bottom, ownD0p>0 = top 
own_tag=0;
ownD0=[0.35;0.3;0.25];ownD0p=[0.6e3;1.2e3;-1.2e3];

% EXPONENTIAL decrease of D0
%--------------------------
exp_tag=0;
D0iniT=0.3;     %Value of D0 near the fault
D0endT=D0top;   %Value of D0 at distD0T
distD0T='0.5*pro_zone'; %give the distance (in m) over which you want the exponential 
                %decrease for the TOP part (!! use string)
D0iniB=0.3;     %Value of D0 near the fault
D0endB=D0bot;   %Value of D0 at distD0B
distD0B='1.0*pro_zone';  %give the distance (in m) over which you want the exponential 
                %decrease for the BOTTOM part (!! use string)
exp_D0

% GAUSSIAN distribution of D0
%---------------------------
gaus_tag=0;
% D0iniT=0.9;     %Value of D0 near the fault
% D0endT=D0top;   %Value of D0 at distD0T
% distD0T='pro_zone'; %give the distance (in m) over which you want the gaussian 
%                 %distribution for the TOP part (!! use string)
% D0iniB=D0bot;     %Value of D0 near the fault
% D0endB=D0bot;   %Value of D0 at distD0B
% distD0B='1e3';  %give the distance (in m) over which you want the gaussian
%                 %distribution for the BOTTOM part (!! use string)
gaus_D0

%% SENSORS POSITIONS

%Give the x and the z coordinates of the synthetic sismometers/sensors
%the xsens values are the actual x coordinates (since variations are along 
%x) but the z coordinates (dzsens) correspond to the distance away from the
%fault
ztp1=[60;240;1020;BOX(4)];
dzsens = [ztp1;ztp1;-ztp1;-ztp1];
xsens = [ones(size(ztp1))*(BOX(1)+1.5e3);ones(size(ztp1))*-3500;...
    ones(size(ztp1))*3500;ones(size(ztp1))*(BOX(2)-1.5e3)];

% field to export for the sensors (you can only give 1)
fieldSen='V';   

%% SIMULATIONS PARAMETERS

TOTtime=4;                  % simulation time duration (in sec)                     
isamp=1;                    % Sampling stride (in number of timesteps) for the sensors.
itd=200;                    % Number of timesteps between snapshots
fieldSnap ='ESVA';          % fields to export in snapshots
ps='F';                     % PostScript (see &SNAP_PS input block)

%Do you want an example with a pure elastic medium? yes=1, no=0;
nodmg_tag=1;

%% 
% =========================================================================
%                              COMPUTATION
% =========================================================================

%% COMPUTE THE PARAMETERS NEEDED FOR THE INPUT FILE

compute_ini_param
% compute_ini_CS_CP

%% BUILT THE MESH

%Create the fault
fault_gen

%ouput names
output_name

% Check if the user want to continue
if checkvarD==1
    disp(' '); disp(['it will erase previous files in ', namefoldD])
	m=input('Do you want to continue? y/n:','s');
	if strcmp(m,'N')==1 || strcmp(m,'n')==1; return; end
end
if checkvarE==1
    disp(' '); disp(['it will erase previous files in ', namefoldE])
	m=input('Do you want to continue? y/n:','s');
	if strcmp(m,'N')==1 || strcmp(m,'n')==1; return; end
end

%Create the mesh
mesh2d_generate

%% OUTPUT FILES & FOLDERS

%export the mesh
mesh2d_write(mesh,[namefoldD,'/',namemesh]);
if nodmg_tag == 1;   mesh2d_write(mesh,[namefoldE,'/',namemesh]); end

%determine dtau for the nucleation prone patch
dtau_comp

%create the output sensors file
write_sensors_file

%create the input file(s) Par.inp
create_inputfile
if nodmg_tag == 1, create_inputfile_nodmg,end


