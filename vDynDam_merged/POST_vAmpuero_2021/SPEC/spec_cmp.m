% 2D Amplification Ratio from SEM2DPACK data
% compared to 1D estimation
% -> contour plots
% H. Rendon and J.P. Ampuero - Sept. 2001
%
% USES:	rock_ampli.m, spec1D.m

clear
disp('This script plots 2D spectra from SEM2DPACK and compares to 1D estimation');
disp('You must modify some parameters in the source file');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOU MUST SET THE FOLLOWING PARAMETERS :

% path to main directory with SEM2DPACK inputs and outputs :
mainpath ='/home/ampueroj/CARACAS/CARACAS_NOV05/Results_2D/';

% list all directories with SEM2DPACK inputs and outputs :
% each directory must contain 
%	Database
%	specf.dat
%	specx.dat
%	xsismos.dat
datapaths = {'pns01/ang146/pwave', 'pns01/ang146/swave', 'pns03/ang22/pwave' };
%datapaths = {'pew/ang37/pwave', 'pew/ang37/swave', 'pew/ang52/pwave', 'pew/ang52/swave' ,...
%	     'pns01/ang146/pwave', 'pns01/ang146/swave', 'pns01/ang22/pwave', 'pns01/ang22/swave' ,...
%	     'pns03/ang146/pwave', 'pns03/ang146/swave', 'pns03/ang22/pwave', 'pns03/ang22/swave' };

% normalization of seismogram amplitude spectra
%	0	no normalization
%	1	by source spectrum
%	2	by rock response
%	3	by a reference simulation
normalize = 2;  

% Scaling options
xscal = 0.001  ; % scaling x-axis
ascal = 1.0   ; % amplitude scaling of 2D spectra

plot1dspec ='N';  % Plot also 1D spectra, only works for two-layer models
% a1/a2 = P velocity (m/s); b1/b2 = S velocity (m/s); ro1/ro2 = density (kg/m^3)
a1 = 2000.; b1 = 900.; ro1 = 2000;               % SEDIMENT
a2 = 3800.; b2 = 2000.; ro2 = 2400;              % ROCK
sedimentfile =[mainpath,'basamento.dat']; 	% file containing sediment bottom (X,Z)

% required if normalize=3 : paths to reference models (ex: without sediments)
refpaths = {'ref1','ref2'}; 

% General plot settings
fmin = 0.5;        %Minimum frequency
fmax = 5.;        %Maximun frequency
R=[0:0.5:3];       %contour scale for 2D
V=[0:0.1:3];       %contour scale for 1D
N    = 40;         %Number of frequencies for 1D spectrum plot

xtitle = 'Location (km)';
ytitle = 'Frequency (Hz)';
tsize = 16; % Font size for titles and axis labels

% END OF USER SETTINGS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datapaths = strcat(mainpath,datapaths);
refpaths = strcat(mainpath,refpaths);
sourcepath = pwd;

wav = 'PS';

disp('Plotting 2D spectra ...')

for ndir = 1:length(datapaths),

  datapath = datapaths{ndir}

  % Read from DataBase :
  % Incident plane wave parameters and title
  [fid,msg] = fopen([datapath,'/DataBase']);
  if fid<0, error([datapath,'/DataBase :' msg]); end
  for n=1:7, plottitle = fgetl(fid); end % move to line #23 -> title
  for n=8:43, line = fgetl(fid); end % move to line #23 -> angle and wave-type
  fclose(fid);
  stuff = sscanf(line,'%f');
  theta = stuff(8); % incidence angle with respect to horizontal X-axis 
  index = stuff(9); % 1 (P-Wave) or 2 (S-Wave)
  plottitle = sprintf('AR %s, %c wave, i=%0.1f',deblank(plottitle),wav(index),theta)

  load([datapath,'/specf.dat']); 
  load([datapath,'/xsismos.dat']); 
  xsismos = xsismos*xscal ;
  load([datapath,'/specx.dat']);

  switch normalize
    case 0
      specx = specx*ascal; 
  
    case 2
      [UXrock,UZrock] = rock_ampli(theta-90,index,a2,b2);
      specx = specx/abs(UXrock) ;
  
    case 3
      specref = load([refpaths{ndir},'/specx.dat']);
      specx = specx ./ specref ;
  end

  kf = find( specf>=fmin & specf<=fmax );
  
   %contourf(xsismos(:,1),specf(kf),specx(:,kf)'); colorbar; 
  [C,H,CF] = contourf(xsismos(:,1),specf(kf),specx(:,kf)',R);
  caxis([R(1) R(end)]) % this forces the color scale to stay fixed !
  colorbar
  xlabel(xtitle,'FontSize',tsize)
  ylabel(ytitle,'FontSize',tsize)
  title(plottitle,'FontSize',tsize)
  
  %---- 1D AR contour plot
  if plot1dspec == 'Y'
      
    hold on;
    disp('Plotting 1D spectra ...');
    depthfile = 'Depthfile';
    spec1D;
    
    if normalize == 2
      UX = UX/UXrock ;
    end
    
    [C,H]=contour(LOC,f,UX,V,'w'); clabel(C,H,V);
    %[C,H]=contour(LOC,f,UX,'w');clabel(C,H)
    
  end
  
  hold off;
  print(gcf,'-depsc2',[datapath,'/spec.eps'])
  
end
