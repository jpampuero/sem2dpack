% *************************************************************************
% This function plot the last snapshost of recorded dynamic damage variables
% such as damage, change in cs or cp, dl/dt, KI, etc... It also plots the slip
% and slip rate evolution along the fault.

% Marion Thomas, last modified August 2018

% CALLS: sem2d_snapshot_read.m; export_fig.m; sem2d_snapshot_plot_km.m

%% ************************************************************************
function [] = plot_snapshots_dyD(datadir,grid,oo,var,T,faultD,faultE,save_tag,namef,opt,writcol,cmap)

%% INPUTS VARIABLES
%datadir:   directory of the data
%grid:      corresponding grid to plot the snapshot
%oo:        potential snapshots
%var:       variables to be plotted
%T:         time for the snapshot
%faultD:	simlations value on the fault for the damage case
%faultE:	simlations value on the fault for the elastic case
%save_tag:  saving option (1=save, 0=just plot, 2=save just .fig)
%namef:     Name of the subfolder for output
%cmap:      colormap for subplot 1
%writcol:	color for the writting

%options main plot
proZ = opt(22);         %process zone size
xb = opt(1:2)./proZ;    %min/max for x-axis
zb = opt(3:4)./proZ;    %min/max for z-axis
as = opt(5);            %aspect ratio between axis (x/z)
mV = opt(6);          %minimum value for the main variable to plot choosen
MV = opt(7);          %maximum value for the main variable to plot choosen
fs = opt(11);           %font size
tplot = opt(23);        %Time at which we want to plot the variable

%option slip and slip rate plots
tstep = opt(8);     %time steps for the slip and slip velocity plots
ym_d = opt(9);      %max value for slip  plot
ym_v = opt(10);     %max value for slip rate plot
as2 = opt(20);      %aspect ratio between axis (x/y) for slip
as3 = opt(21);      %aspect ratio between axis (x/y) for slip rate

%comments plot
optplot=[opt(1:4)./proZ,opt(11:19)];

%% PLOT PARAMETERS

%find snpashot to plot
o=find(abs(T-tplot)== min(abs(T-tplot)));
touput=num2str(round(10*T(oo(o)+1))/10);

%Launch figure
v_handle =figure('Position',[200 500 1600 750],'PaperOrientation','landscape');  
f=v_handle;
set(f,'PaperOrientation','landscape');
set(f,'PaperUnits','normalized');
set(f,'PaperPosition', [0 0 1 1]);
namefile=[namef,var,'_slip_sliprate'];

%% DOWNLOAD

%downlaod
dyD = sem2d_snapshot_read('dyD',oo(o),datadir);

%Find the variable to plot
if strcmp('D',var)==1
    v=dyD.D;
    vn='Damage'; vnC='D';
elseif strcmp('cs',var) == 1
    v=dyD.cs;
    cs0=max(max(max(v)));
    v=100*(cs0-v)/cs0;
    vn='Cs variation'; vnC='%';
elseif strcmp('cp',var) == 1
    v=dyD.cp;
    cp0=max(max(max(v)));
    v=100*(cp0-v)/cp0; 
    vn='Cp variation';vnC='%';
elseif strcmp('dldt',var) == 1
    v=dyD.dldt;
    vn='dldt';vnC='dldt';
elseif strcmp('R',var) == 1
    v=dyD.R;
    vn='Regime';vnC='R';
else
    'not a dyD variables'
end

%comments plot
disp(' ')
disp([vn ' is plotted with cumulative slip and slip rate'])

%% CORRECTION FOR CP CS

%Find the variable to plot
D=dyD.D;
idD=find(D<=max(opt(12),opt(16)));

if (strcmp('cs',var) == 1)
    v(idD)= 0;
elseif strcmp('cp',var) == 1
    v(idD)= 0;
end

%% PLOT SNAPSHOT

ax=subplot(2,2,[1,2]); hold on
if (save_tag==0),sem2d_snapshot_plot_norm(v,grid,proZ,[mV MV]);end
plot(faultD.x./proZ,faultD.z./proZ,'-w','Linewidth',0.5)
colormap(cmap);
t=colorbar('peer',gca);
set(get(t,'xlabel'),'String', vnC,'Fontweight','bold')
annotation_plot('',touput,optplot,writcol,as,ax,2,2);

%% Slip and slip velocity

%Color
ctab=[0.5 0 0.5;1 0 0;1 0.5 0;1 0.8 0;0 0.9 0;0 0.4 0.4;0 0 0.4;0 0 1;0 0.9 0.9];

%Options
maxtime=tplot;
stp=round(size(faultD.t,2)*tstep/max(faultD.t));
stpE=round(size(faultE.t,2)*tstep/max(faultE.t));
id=find(abs(faultD.t-maxtime)== min(abs(faultD.t-maxtime)));
idE=find(abs(faultE.t-maxtime)== min(abs(faultE.t-maxtime)));

%LEFT slip
subplot(2,2,3)
hold on; for i=1:stp:id; plot(faultD.x./proZ,faultD.d(:,i),'color',[0 0.7 1],'Linewidth',2);end
hold on; for i=1:stpE:idE; plot(faultE.x./proZ,faultE.d(:,i),'k','Linewidth',1);end
text(xb(1)+((xb(2)-xb(1))*0.38),ym_d-ym_d*0.05,[' time step = ', num2str(round(100*faultD.t(1+stp))/100),' s'],'FontSize', fs);
xlabel('Distance/R_0','FontSize', fs)
ylim([0 ym_d]); ylabel('Cumulative slip (m)','FontSize', fs)
xlim(xb);
ax = gca;
ax.FontSize = fs;
daspect([1 as2 1])

%RIGHT slip rate
subplot(2,2,4)
hold on; for i=1:stpE:idE; plot(faultE.x./proZ,faultE.v(:,i),'k','Linewidth',1);end
hold on; for i=1:stp:id; plot(faultD.x./proZ,faultD.v(:,i),'-k','MarkerSize',3);ylim([0 ym_v]);end
for j=1:size(ctab,1)
    plot(faultD.x./proZ,faultD.v(:,i),'o-','MarkerSize',3,'color',ctab(j,:));
    i=i-stp;
end
text(xb(1)+((xb(2)-xb(1))*0.38),ym_v-ym_v*0.05,[' time step = ', num2str(round(100*faultD.t(1+stp))/100),' s'],'FontSize', fs);
xlabel('Distance/R_0','FontSize', fs)
ylabel('Slip rate (m/s)','FontSize', fs)
xlim(xb);
ax = gca;
ax.FontSize = fs;
daspect([1 as3 1])


%% SAVE FIGURES

%for illustrator
if (save_tag>=1) 
    
    %Save in eps format
    if (save_tag==1),export_fig(namefile,v_handle,'-eps'),end

    %Plot the screenshot
    ax=subplot(2,2,1:2); hold on
    sem2d_snapshot_plot_norm(v,grid,proZ,[mV MV]);
    plot(faultD.x./proZ,faultD.z./proZ,'-w','Linewidth',0.5)
    t=colorbar('peer',gca);
    set(get(t,'xlabel'),'String', 'D','Fontweight','bold')
    xlim(xb); ylim(zb)
    daspect([1 as 1])

	%Save in png format
    if (save_tag==1), export_fig(namefile,gcf,'-png','-r340');    end

    %Comments
    annotation_plot('',touput,optplot,writcol,as,ax,2,2);

	%Save in matlab format
    saveas(v_handle,[namefile,'.fig'],'fig')
end

%% Slip and slip velocity2

v_handle2 =figure('Position',[200 500 350 140],'PaperOrientation','landscape');  

%LEFT slip
subplot(1,2,1)
hold on; for i=1:stp:id; plot(faultD.x./proZ,faultD.d(:,i),'color',[0 0.7 1],'Linewidth',2);end
hold on; for i=1:stpE:idE; plot(faultE.x./proZ,faultE.d(:,i),'k','Linewidth',1);end
ylim([0 ym_d/4]); %ylabel('Cumulative slip (m)','FontSize', fs)
xlim([xb(1)+(xb(2)-xb(1))/24 xb(1)+(xb(2)-xb(1))/6]);
ax = gca; ax.FontSize = fs;
daspect([1 as2 1])
subplot(1,2,2)
hold on; for i=1:stp:id; plot(faultD.x./proZ,faultD.d(:,i),'color',[0 0.7 1],'Linewidth',2);end
hold on; for i=1:stpE:idE; plot(faultE.x./proZ,faultE.d(:,i),'k','Linewidth',1);end
ylim([0 ym_d/4]); %ylabel('Cumulative slip (m)','FontSize', fs)
xlim([xb(2)-(xb(2)-xb(1))/6 xb(2)-(xb(2)-xb(1))/24]);
ax = gca;ax.FontSize = fs;
daspect([1 as2 1])

%for illustrator
if (save_tag>=1) 
    namefile=[namef,vn,'_slip_sliprate_zoom'];
    if (save_tag==1),export_fig(namefile,v_handle2,'-eps'),end
end
