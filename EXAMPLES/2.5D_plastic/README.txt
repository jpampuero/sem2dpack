---------------------------------------------------------------------------------------
EB_initial_stresses_and_cohesion.m is used to generate initial values for sem2Dpack input file Par.inp
---------------------------------------------------------------------------------------

E.g. in psi_45_S_0.56_CF_0.63_W_10/Par.inp:

#---- Material parameters --------------
&MATERIAL tag=1, kind='PLAST'  /
&MAT_PLASTIC cp=5770.d0, cs=3330.d0, rho=2705.d0, phi = 30.d0, coh = 9.669501d6, Tv = 0.0561d0,
                 e0 = -4.162378e-04, -4.162378e-04, 3.504802e-04 /

#----- Boundary conditions ---------------------
&BC_DEF  tags = 5,6 , kind = 'DYNFLT' /
&BC_DYNFLT friction='SWF','TWF', Tn=-50d6, Tt=2.102564d7 /
&BC_DYNFLT_SWF Dc=2d0, MuS=0.6d0, MuD=0.1d0 /
&BC_DYNFLT_TWF kind=1, MuS=0.6d0, MuD=0.1d0, Mu0=0.6d0,
               X=0.d0, Z=0.d0, V=0.333d3, L=0.1665d3, T=60d0 /
			   

where coh, e0, and Tt values are obtained by running EB_initial_stresses_and_cohesion.m
with:

%% Input values
psi    = 45;              % direction of max principal stress to fault (degrees)
S      = 0.56;            % strength excess parameter = (taus-tau0)/(tau0-taud)
szz    = -50e6;           % fault normal stress (sigma)
W      = 10e3;
CF     = 0.63;            % closeness to failure (Mohr circle radius / distance to yield surface)

%% Material properties
cp     = 5770;            % p-wave velocity
cs     = 3330;            % s-wave velocity
rho    = 2705;            % density
phi    = 30;              % internal friction angle (degrees)

%% Frictional properties (SW friction)
mus    = 0.6;            % static friction coefficient
mud    = 0.1;            % dynamic friction coefficient
Dc     = 2;              % slip-weakening distance

coh - Cohesion, e0 - Initial strain (exx, ezz, exz), Tt - Fault shear stress

---------------------------------------------------------------------------------------
EB_plot_fault_and_energies.m is used to plot the simulation results after running sem2dsolve Par.inp
---------------------------------------------------------------------------------------