% SPEC2D_CPLOT Distance-frequency contour plot of 2D Amplification Ratio (AR)
%
% AR is the spectral ratio between SEM2DPACK seismograms with plane wave source
% and the reference seismogram at the surface of a uniform half-space with
% vertically incident plane wave source and same source time function as in SEM2DPACK.
%
% The valid range of frequencies [f1,f2] for the plot depends on the source time
% function of the SEM2DPACK simulation: it must contain significant source 
% spectral power, say more than 1% of its spectral peak. For a Ricker source
% with dominant frequency f0, [f1,f2] is typically [0.1,2.5]*f0.
%	
% AUTHORS	H. Rendon and J.P. Ampuero - Sept. 2001 
% LAST UPDATE 	Thu Aug 14 17:09:51 PDT 2008 (Ampuero)

% TO DO		+ define the reference seismogram with same oblique incidence as in SEM2DPACK
%		  (modify ascal, using rock_ampli.m)
%		+ improve the spectral amplitude estimation, needs at least some tapering

plottitle = 'AR at Los Palos Grandes, X component, P wave, i=60 N';
datapath = '.';	% directory containing SEM2DPACK output files
comp = 'y';	% select seismogram component (x, y or z)
f1 = 0.1  ; 	% minimum frequency, typically 0.1*f0 for a Ricker source
f2 = 2.5; 	% maximum frequency, typically 2.5*f0 for a Ricker source
nf = 40 ; 	% minimum number of frequencies in the range [f1,f2]
ascal = 2;	% surface amplification factor for vertical incidence on a uniform half-space
xscal = 1000;	% scaling x-axis
xtitle = 'Location (km)';
ytitle = 'Frequency (Hz)';
tsize = 14; % Font size for titles and axis labels

data = sem2d_read_seis(datapath);
u = getfield(data,['u' lower(comp)]);
src = load([datapath '/SourcesTime_sem2d.tab']);

nfft = max( data.nt, ceil(nf/data.dt/(f2-f1)));
nfft = 2^nextpow2(nfft);
u = fft(u,nfft);
src = fft(src(:,2),nfft);

freq = [0:nfft/2]/(nfft*data.dt);
n1 = find(freq>=f1, 1,'first');
n2 = find(freq<=f2, 1,'last');
ar = abs(u(n1:n2,:))./repmat(abs(src(n1:n2)),1,size(u,2))/ascal; 

contourf(data.x/xscal, freq(n1:n2), ar, 30,'LineStyle','none'); 
ca=caxis;
caxis([0 ca(2)]);
colorbar; 
xlabel(xtitle,'FontSize',tsize);
ylabel(ytitle,'FontSize',tsize);
title(plottitle,'FontSize',tsize);
