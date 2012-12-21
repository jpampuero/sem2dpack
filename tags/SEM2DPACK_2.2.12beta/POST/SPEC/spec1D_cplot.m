% 1D single sediment layer Amplification Ratio
% Contour plot
% H.Rendon Set. 2001

clear;
disp('1D amplification spectral estimation');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define Model Parameters
a1 = 1775.; b1 = 1025.; ro1 = 1800;              %BASEMENT Half Space
a2 = 4000.; b2 = 2300.; ro2 = 2200;              %SEDIMENTARY Layer

% Plot settings
fmin = 0.5;       %Minimum frequency
fmax = 5.;        %Maximun frequency
N    = 50;       %Number of frequencies
xscal = 0.001 ;
% Look also at shifting of x-axis, some lines down
% and fix contour options (at the end of this file)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Define the diferent depths along the Valley
depthfile = input('Enter depth file : ','s');

% Define the incident wave
theta = input('Enter incidence angle with respect X-axis => ');      %Angle of incidence
index = input('(1) P-Wave,    (2) S-Wave   ');

% Compute
spec1D ;

% Contour plot
V=[0 : 0.5 : 4]; [C,Hc]=contour(LOC,f,UX,V,'k');clabel(C,Hc,V);
 %[C,H]=contour(LOC,f,UX,'w');clabel(C,H);
