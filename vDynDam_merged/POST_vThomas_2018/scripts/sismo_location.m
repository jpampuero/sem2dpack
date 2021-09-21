% ************************************************************************
% This function plot the particules velocity for the different "seismometers"
% (point outputs, definied in the input file)

% Marion Thomas, last modified october 2018

% CALLS: sem2d_snapshot_read.m; export_fig.m; sem2d_snapshot_plot_norm.m

%% ************************************************************************
function [] = sismo_location(datadir,grid,oo,T,Pcoord,data,fault,tend,save_tag,namef,opt,writcol,cmap)

%% INPUTS VARIABLES

%datadir:	directory of the data
%grid:      corresponding grid to plot the snapshot
%oo:        potential snapshots
%T:         time
%Pcoord:	seismographs locations
%data:      seismograms (potentially including several simulations)
%save_tag:  saving option (1=save, 0=just plot, 2=save just .fig)
%namef      field to be plot
%field      definie which component (z or x) should be plotted
%cmap:      colormap for subplot 1
%writcol:	color for the writting

%options
proZ = opt(15);         %process zone size
xb = opt(1:2)./proZ;    %min/max for x-axis
zb = opt(3:4)./proZ;    %min/max for z-axis
as = opt(5);            %min/max for x-axis, y-axis and the aspect ratio between axis (x/y)

%comments plot
optplot=[opt(1:4)./proZ,opt(6:14)];

%% PLOT PARAMETERS I

%find snpashot to plot
% o=find(abs(T-tend(1))== min(abs(T-tend(1))));
tplottrue=T(end);
touput=num2str(round(10*tplottrue)/10);

%selected seismograms
nsS=size(Pcoord,1);         %number of selected seismograms

%curve color
ctab=[0 0 0;1 0 0; 0 0 1; 0 1 0; 1 0 1; 1 1 0; 0 1 1 ; 0.5 0.5 0.5];

%% DOWNLOAD DATA

%snapshot
dyD = sem2d_snapshot_read('dyD',oo(end),datadir);
v=dyD.D;

%number of simulations
Nsimu=size(data,1);
s(1).nt=0;
s(1).ns=0;

for i=1:Nsimu
%x and z coordinates
    s(i).x=data(i).x;           %x-coord of seismograms
    s(i).z=data(i).z;           %z-coord of seismograms
    
%fault coordinates    
    f(i).x=fault(i).x;           %x-coord of fault
    f(i).z=fault(i).z;           %z-coord of fault

%% SELECTED SEISMOGRAMS

    %find the unique x coordinates
    xsens=Pcoord(:,1);
    xloc1=unique(Pcoord(:,1));
    %only keep the ones that correspond to a x coordinates one the fault
    xloc=xloc1(find(xloc1>=f(i).x(1) & xloc1 <= f(i).x(end)));
    %find the corresponding z coordinates on the fault
    z_xloc=zeros(size(xloc));
    for j=1:numel(xloc)
        dxloc=abs(xloc(j)-f(i).x);
        z_xloc(j)=f(i).z(find(dxloc==min(dxloc)));
    end
    %Define the z coordinates of the sismometers
    zsens=Pcoord(:,2);
    for j=1:numel(xloc)
        idS=find(xsens == xloc(j));
        zsens(idS)=zsens(idS)+z_xloc(j);
    end

    %Find the seismograms the closest to the location given by the user
    for k=1:nsS
        matfind_p=abs(s(i).x-xsens(k))+abs(s(i).z-zsens(k));
        s(i).ind(k)= find(min(matfind_p) == matfind_p);  %n0 of the selected seismograms
    end
end

%% SEISMO LOCATION

%Launch figure
v_handle2=figure('Position',[0 500 1600 1000],'PaperOrientation','landscape');hold on
h=v_handle2;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
namefile2=[namef,'seismo_loc'];

%Damage
ax=subplot(2,1,1); hold on
if (save_tag==0),sem2d_snapshot_plot_norm(v,grid,proZ,[0 1]);end
colormap(cmap);
t=colorbar('peer',gca);
set(get(t,'xlabel'),'String', 'D','Fontweight','bold')
xlim(xb); ylim(zb)
daspect([1 as 1])
annotation_plot('',touput,optplot,writcol,as,ax,1,1)

%Points
for i=1:Nsimu
	plot(s(i).x(s(i).ind)./proZ, s(i).z(s(i).ind)./proZ,'o',...
	'MarkerFaceColor',ctab(i,:),'MarkeredgeColor','w')
end
xlim(xb); ylim(zb)

%for illustrator
if (save_tag>=1) 
    
	%Save in eps format
    if (save_tag==1),export_fig(namefile2,v_handle2,'-eps'),end

    %Plot the screenshot
    ax=subplot(2,1,1); hold on
    sem2d_snapshot_plot_norm(v,grid,proZ,[0 1])
    t=colorbar('peer',gca);
    set(get(t,'xlabel'),'String', 'D','Fontweight','bold')
    xlim(xb); ylim(zb)
    daspect([1 as 1])

	%Save in png format
    if (save_tag==1), export_fig(namefile2,gcf,'-png','-r340');    end

    %Comments
    annotation_plot('',touput,optplot,writcol,as,ax,1,1)
    
    %Points
    for i=1:Nsimu
        plot(s(i).x(s(i).ind)./proZ, s(i).z(s(i).ind)./proZ,'o',...
        'MarkerFaceColor',ctab(i,:),'MarkeredgeColor','w')
    end

	%Save in matlab format
    saveas(v_handle2,[namefile2,'.fig'],'fig')
end

