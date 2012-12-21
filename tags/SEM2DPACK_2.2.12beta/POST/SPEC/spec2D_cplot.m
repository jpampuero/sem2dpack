% 2D Amplification Ratio from SPECFEM data
% H. Rendon and J.P. Ampuero - Set. 2001 

xscal = 0.001  ;% scaling x-axis
ascal = 0.5     ;% amplitude scaling
xtitle = 'Location (km)';
ytitle = 'Frequency (Hz)';
plottitle = 'AR at Los Palos Grandes, X component, P wave, i=60 N';
tsize = 18; % Font size for titles and axis labels

disp('Plots spectra map');
datapath = input('Enter data path : ','s'); sourcepath = pwd;
cd(datapath); load specf.dat; load xsismos.dat; load specx.dat; cd(sourcepath);
specx = specx*ascal; xsismos = xsismos*xscal ;
nf = size(specf,1);

f1 = 5  ; % first frequency to be plotted
f2 = nf ; % last frequency to be plotted

contourf(xsismos(:,1),specf(f1:f2),specx(:,f1:f2)'); colorbar; 
xlabel(xtitle,'FontSize',tsize);ylabel(ytitle,'FontSize',tsize);
title(plottitle,'FontSize',tsize);
