% Script to download data from the output folder

% Marion Thomas last modified July 2018

%CALLS: 

%==========================================================================
%% 
%load the grids
gridD=sem2d_read_specgrid(datadirD);
gridE=sem2d_read_specgrid(datadirE);
gridF=sem2d_read_specgrid(datadirF);

%Fault
faultD = sem2d_read_fault([datadirD,'Flt05']);
faultD.t =0:faultD.dt:(faultD.nt-1)*faultD.dt;
% timeD=linspace(0,(faultD.nt-1)*faultD.dt,faultD.nt-1);
faultE = sem2d_read_fault([datadirE,'Flt05']);
faultE.t =0:faultE.dt:(faultE.nt-1)*faultE.dt;
% timeE=linspace(0,(faultE.nt-1)*faultE.dt,faultE.nt-1);
faultF = sem2d_read_fault([datadirF,'Flt05']);
faultF.t =0:faultF.dt:(faultF.nt-1)*faultF.dt;
% timeF=linspace(0,(faultF.nt-1)*faultF.dt,faultF.nt-1);
   
%Sismometers
%locations has been definied in the input 'sensors.ini'
sismo = sem2d_read_seis(datadirD); %damage + rough fault
timesis=linspace(0,(sismo.nt-1)*sismo.dt,sismo.nt);
sismoE = sem2d_read_seis(datadirE); %elasticity + rough fault
timesisE=linspace(0,(sismoE.nt-1)*sismoE.dt,sismoE.nt);
sismoF = sem2d_read_seis(datadirF);%elasticity + flat fault
timesisF=linspace(0,(sismoF.nt-1)*sismoF.dt,sismoF.nt);

%compute the value of the process zone
cs=[infoD(1,2),infoD(end,2)];
rho=[infoD(1,4),infoD(end,4)];
nu=[infoD(1,5),infoD(end,5)];
mu=rho.*cs.*cs;
Szz=max(faultD.sn0);
if (ndof==2);mu_s=mu./(1-nu);else mu_s=mu;end
pro_zoneTmp = (9*pi/32)*mu_s/((-MuS+MuD)*Szz);
pro_zone = max(pro_zoneTmp);
disp(['process zone:         ',num2str(pro_zone/1e3), ' km'])

%geometry info
zuni=unique(round(gridD.coorg(:,2)));
Mz=max(gridD.coorg(:,2));
mz=min(gridD.coorg(:,2));
xuni=unique(gridD.coorg(:,1));
Mx=max(gridD.coorg(:,1));
mx=min(gridD.coorg(:,1));

%Resolution
hD=round((Mx-mx)/(numel(xuni)-1));
% hD=(faultD.x(end)-faultD.x(1))/(size(faultD.x,1)-1);
hE=(faultE.x(end)-faultE.x(1))/(size(faultE.x,1)-1);
hF=(faultD.x(end)-faultF.x(1))/(size(faultF.x,1)-1);

%Simulation duration:


disp(' ');
disp('========================================');
disp('Information about the simulation')
disp('========================================');
disp(['Domain size: ',num2str(Mx-mx), 'x',num2str(Mz-mz),' m'])
disp(['x-coordinates: [',num2str(mx), ',',num2str(Mx),'] in m'])
disp(['z-coordinates: [',num2str(mz), ',',num2str(Mz),'] in m'])
disp(['resolution: ' num2str(hD),' m'])
disp(['Total simulation duration = ', num2str(simuD),'s']);
