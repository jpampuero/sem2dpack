% 2D Amplification Ratio from SEM2DPACK data
% compared to 1D estimation
% -> contour plots
% H. Rendon and J.P. Ampuero - Sept. 2001

home; clear;
disp('This script plots 2D spectra from SEM2DPACK and compares to 1D estimation');
disp('You must modify some parameters in the source file');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOU MUST SET THE FOLLOWING PARAMETERS :

% Data Paths
% datapath = '/home/rocabado/PalosGrandes/Mod-P/P-60N/'  ;          % -> results from SEM2DPACK
% sedimentfile ='/home/rocabado/PalosGrandes/Mod-P/PNS4-Depth.txt';   % -> sediment (X,Z)

datapath = '/home/rocabado/PNS02/Mod-P/P-60N/'  ;          % -> results from SEM2DPACK
sedimentfile ='/home/rocabado/PNS02/Mod-P/ns02d.dat';     % -> sediment (X,Z)

% Normalize spectra by rock response or not ?  ('Y' or 'N')
rocknormalize = 'N';

% Normaliser par une autre simu:
simunormalize = 'Y';
datapath2 = ' '   

% Incident plane wave parameters
theta = 150; % incidence angle with respect to horizontal X-axis 
index = 1 ;    % (1) P-Wave,    (2) S-Wave   

% Model Parameters
a1 = 2000.; b1 = 900.; ro1 = 1800;               % SEDIMENT
a2 = 3800.; b2 = 2000.; ro2 = 2300;              % ROCK

% Scaling options
xscal = 0.001  ; % scaling x-axis
ascal = 0.5   ; % amplitude scaling of 2D spectra

% General plot settings
plot1D = 'N'
fmin = 0.5;        %Minimum frequency
fmax = 2.5;        %Maximun frequency
N    = 40;         %Number of frequencies for 1D spectrum plot
R=[0:0.4:6];       %contour scale for 2D
V=[1 : 1 :10];     %contour scale for 1D

xtitle = 'Location (km)';
ytitle = 'Frequency (Hz)';
plottitle = 'AR at San Bernardino, P wave, i=-60 N';
tsize = 16; % Font size for titles and axis labels

% END OF USER SETTINGS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(plottitle)
disp('Plotting 2D spectra ...')
sourcepath = pwd;

cd(datapath); 
load specf.dat; load xsismos.dat; load specx.dat;

if simunormalize == 'Y'
  cd(datapath2); 
  specref = load('specx.dat');
  specx = specx/specref ;
end
  
cd(sourcepath);
specx = specx*ascal; xsismos = xsismos*xscal ;
nf = size(specf,1);

if rocknormalize == 'Y'
  rock_ampli;
  specx = specx/UXrock ;
end

f1 = max(find( specf<fmin ))+1 ;
f2 = min(find( specf>fmax ))-1 ;

 %contourf(xsismos(:,1),specf(f1:f2),specx(:,f1:f2)'); colorbar; 
[C,H,CF] = contourf(xsismos(:,1),specf(f1:f2),specx(:,f1:f2)',R);colorbar;
xlabel(xtitle,'FontSize',tsize);ylabel(ytitle,'FontSize',tsize);
title(plottitle,'FontSize',tsize);

if plot1D == 'Y'

hold on;
disp('Plotting 1D spectra ...');
depthfile = 'Depthfile';
spec1D;

if rocknormalize == 'Y'
  UX = UX/UXrock ;
end

%---- 1D AR contour plot
[C,H]=contour(LOC,f,UX,V,'w'); clabel(C,H,V);
%[C,H]=contour(LOC,f,UX,'w');clabel(C,H)

hold off;

end
