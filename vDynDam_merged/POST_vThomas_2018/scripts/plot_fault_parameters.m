% ************************************************************************
% This function plot the last snapshost of recorded dynamic damage variables
% such as damage, change in cs or cp, dl/dt, KI, etc... It also plots the slip
% and slip rate evolution along the fault.

% Marion Thomas, last modified May 2017

% CALLS: sem2d_snapshot_read.m; export_fig.m; sem2d_snapshot_plot_km.m

%% ************************************************************************
function [] = plot_fault_parameters(data,dataE,save_tag,namef,opt)

% data=faultD;
% dataE=faultE;
% namef=namefoldD;
% opt=options;
% 
% close all; clc;

%% INPUTS VARIABLES
%datadir:   directory of the data
%grid:      corresponding grid to plot the snapshot
%oo:        potential snapshots
%varN:      variable to be plotted
%T:         time for snapshots
%data:      simlations value on the fault for the damage case
%dataE:     simlations value on the fault for the elastic case
%save_tag:  saving option (1=save, 0=just plot, 2=save just .fig)
%namef:     folder
%cmap:      colormap for subplot 1

%options
proZ = opt(11);         %process zone size
xb = opt(1:2)./proZ;    %min/max for x-axis
ndt = opt(3);           %number of iscochrone
ym_d = opt(4);          %max slip for slip  plot
ym_v = opt(5);          %max slip for slip rate plot
ym_sn = opt(6:7);       %max slip for slip rate plot
ym_st = opt(8:9);       %max slip for slip rate plot
fs = opt(10);           % font size
tstart = opt(12);       %Time at which we want to start plotting the variable
tend = opt(13);         %Time at which we want to end plotting the variable

%initial value
sn0=abs(min(data.sn0));
st0=min(data.st0);
snE0=abs(min(dataE.sn0));
stE0=min(dataE.st0);

%% PLOT PARAMETERS

%Launch figure
v_handle =figure('Position',[200 500 1600 750],'PaperOrientation','landscape');  
f=v_handle;
set(f,'PaperOrientation','landscape');
set(f,'PaperUnits','normalized');
set(f,'PaperPosition', [0 0 1 1]);
namefile=[namef,'faultparam'];

% Slip and slip velocity color code
ctab=[0.5 0 0.5;1 0 0;1 0.5 0;1 0.8 0;0 0.9 0;0 0.4 0.4;0 0 0.4;0 0 1;0 0.9 0.9];

%comment
disp(' ')
disp(['display fault parameters over time'])

%% Slip and slip velocity

%Options
id=find(abs(data.t-tend)== min(abs(data.t-tend)));
idE=find(abs(dataE.t-tend)== min(abs(dataE.t-tend)));
id0=find(abs(data.t-tstart)== min(abs(data.t-tstart)));
idE0=find(abs(dataE.t-tstart)== min(abs(dataE.t-tstart)));
stp=floor(size(data.t(id0:id),2)./ndt);
stpE=floor(size(dataE.t(idE0:idE),2)./ndt);

%slip
% ax=subplot(4,1,1);
ax=subplot(2,1,1);
hold on; 
for i=idE0:stpE:idE; plot(dataE.x./proZ,dataE.d(:,i),'k','Linewidth',1);...
        %display(['tE=',num2str(timeE(i))]);
end
j=1;
for i=id0:stp:id
    plot(data.x./proZ,data.d(:,i),'-','MarkerSize',2,'color',ctab(j,:),'Linewidth',1);j=j+1;
%     display(['t=',num2str(time(i))]);
end
text(xb(1)+((xb(2)-xb(1))*0.38),ym_d-ym_d*0.05,[' time step = ', num2str(round(100*data.t(1+stp))/100),' s'],'FontSize', fs);
set(ax,'Xticklabel',[]);
xlim(xb); ylim([0 ym_d]); 
ylabel('Cumulative slip (m)','FontSize', fs)
ax = gca;
ax.FontSize = fs;

%slip rate
% ax=subplot(4,1,2);
ax=subplot(2,1,2);
hold on; for i=id0:stpE:idE; plot(dataE.x./proZ,dataE.v(:,i),'k','Linewidth',1);end
j=1;
for i=id0:stp:id
    plot(data.x./proZ,data.v(:,i),'-','MarkerSize',2,'color',ctab(j,:),'Linewidth',1);
    j=j+1;
end
% set(ax,'Xticklabel',[]);
xlabel('Distance/R_0','FontSize', fs)
ylabel('Slip rate (m/s)','FontSize', fs)
xlim(xb); ylim([0 ym_v]); 
ax = gca;
ax.FontSize = fs;

% %normal stress
% ax=subplot(4,1,3);
% hold on; for i=idE0:stpE:idE; plot(dataE.x./proZ,dataE.sn(:,i)./snE0,'k','Linewidth',1);end
% j=1;
% for i=id0:stp:id
%     plot(data.x./proZ,data.sn(:,i)./sn0,'-','MarkerSize',2,'color',ctab(j,:),'Linewidth',1);
%     j=j+1;
% end
% set(ax,'Xticklabel',[]);
% xlim(xb);ylim([ym_sn]); 
% ylabel('\sigma_{zz}(t)/\sigma_{zz}(0)','FontSize', fs)
% ax = gca;
% ax.FontSize = fs;
% 
% %shear stress
% ax=subplot(4,1,4);
% hold on; 
% for i=idE0:stpE:idE; plot(dataE.x./proZ,dataE.st(:,i)./stE0,'k','Linewidth',1);end
% j=1;
% for i=id0:stp:id
%     plot(data.x./proZ,data.st(:,i)./st0,'-','MarkerSize',2,'color',ctab(j,:),'Linewidth',1);
%     j=j+1;
% end
% xlabel('Distance/R_0','FontSize', fs)
% xlim(xb); ylim([ym_st]); 
% ylabel('\sigma_{xz}(y)/\sigma_{xz}(0)','FontSize', fs)
% ax = gca;
% ax.FontSize = fs;


%% SAVE FIGURES

%for illustrator
if (save_tag>=1) 
    saveas(v_handle,[namefile,'.fig'],'fig')
    if (save_tag==1),export_fig(namefile,v_handle,'-eps'),end;
end


