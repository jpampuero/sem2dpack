% ************************************************************************
% This function plot the last snapshost of recorded dynamic damage variables
% such as damage, change in cs or cp, dl/dt, KI, etc... It also plots the slip
% and slip rate evolution along the fault.

% Marion Thomas, last modified May 2017

% CALLS: sem2d_snapshot_read.m; export_fig.m; sem2d_snapshot_plot_km.m
% clc;close all
%% ************************************************************************
function [] = plot_profil_dyD(datadir,grid,oo,T,fault,save_tag,namef,opt,writcol,cmap)
% datadir=datadirD;
% grid=gridD;
% oo=ooD;
% T=TD;
% fault=faultD;
% namef=namefoldD;
% opt=options;
% writcol='white';
% cmap=cmapD;

%% INPUTS VARIABLES

%datadir:   directory of the data
%grid:      corresponding grid to plot the snapshot
%oo:        potential snapshots
%var:       variable to be plotted
%T:         time
%fault:     simlations value on the fault for the damage case
%name:      Reference name to the simulation
%save_tag:  saving option (1=save, 0=just plot, 2=save just .fig)
%namef:     folder
%cmap:      colormap for subplot 1
%writcol:	color for the writting

%options snapshot plot
proZ = opt(18);         %process zone size
xb = opt(1:2)./proZ;    %min/max for x-axis
zb = opt(3:4)./proZ;    %min/max for z-axis
as = opt(5);            %aspect ratio between axis (x/z)
mV = opt(6);            %minimum value for the main variable
MV = opt(7);            %maximum value for the main variable
fs = opt(8);            %font size
np=opt(17);             %number of profiles
tplot = opt(19);        %Time at which we want to plot the variable
smooth_op = opt(20);    %smoothing option for snapshot: 1 for yes 0 for no.

%options profiles plot
fac=1/np*10;
D01 = round(10*opt(9))/10*fac;    %initial damage for lower part

%comments plot
optplot=[opt(1:4)./proZ,opt(8:16)];
disp(' ')
disp(['Damage snapshot at ', num2str(tplot), 's is plotted']) 

%% PARAMETERS

%grid resolution
h=(fault.x(end)-fault.x(1))/(size(fault.x,1)-1);

%find snpashot to plot
o=find(abs(T-tplot)== min(abs(T-tplot)));
touput=num2str(round(10*T(oo(o)+1))/10);

%Launch figure
v_handle=figure('Position',[200 500 1600 750],'PaperOrientation','landscape');hold on
f=v_handle;
objf=findobj('type','figure');numfig=length(objf);
set(f,'PaperOrientation','landscape');
set(f,'PaperUnits','normalized');
set(f,'PaperPosition', [0 0 1 1]);

%name files and folder
namef=[namef,'/Damage_profiles/'];
if exist(namef) == 7; else; mkdir(namef); end
namefile=[namef,'Damage_profiles_snap',num2str(oo(o))];


%% DOWNLOAD

%downlaod
dyD = sem2d_snapshot_read('dyD',oo(o),datadir);
dyD0 = sem2d_snapshot_read('dyD',oo(1),datadir);

%Find the variable to plot
v=dyD.D;
v0=dyD0.D;

%% FIND THE PROFILS

% Define the location of the profils

%Distance between profiles
edge=(0.9e3)/proZ;
dx=abs((xb(2)-2*edge-xb(1)))/(np-1);
xp=xb(1)+edge:dx:xb(2)-edge;
dxp=mean(xp(2:end)-xp(1:end-1));  
disp(['with profiles across the fault (every ',num2str(round(dxp*10)/10), ' R_0)']);

%Matrix rearrangement
data_matrix_sorting;

%For profiles
c0=reshape(v0(wrk.indx(:,:)),ll,cc)';
x2=mean(x)'/proZ;
z2=mean(z)'/proZ;
c2=mean(field)';
c02=mean(c0)';

%declare variables
Ny2=round(Nz*(4/(Nz*Nx/grid.nelem)));
Dpro=zeros(Ny2,np);
D0pro=zeros(Ny2,np);
Zpro=zeros(Ny2,np);
Dstdpro=zeros(Ny2,np);

% Find profiles 
for j=1:np
    minloc=abs(x2-xp(j));
    int=find(minloc ==  min(minloc));
%     du=numel(int)/numel(unique(y2(int)));
    du=numel(int)/Ny2;
    intDstd=find(x2<=xp(j)+dxp/2 & x2>=xp(j)-dxp/2);
    c2slc=c2(intDstd);
    Dslc=reshape(c2slc,numel(intDstd)/(Ny2),Ny2);
    Dstdpro(:,j)=smooth(std(Dslc)');
    if du==2
        id1=sort([1:du^2:numel(int),2:du^2:numel(int)]);
        id2=sort([du+1:du^2:numel(int),du+2:du^2:numel(int)]);
        Dpro(:,j)=(c2(int(id1))+c2(int(id2)))/2;
        D0pro(:,j)=(c02(int(id1))+c02(int(id2)))/2;
        Zpro(:,j)=(z2(int(id1))+z2(int(id2)))/2;
    elseif du ==1
        Dpro(:,j)=c2(int);
        D0pro(:,j)=c02(int);
        Zpro(:,j)=z2(int);
    end
end

%normalization
Dprofil=Dpro.*fac;
D0profil=D0pro.*fac;
Dstd=Dstdpro.*fac;
Zprofil=Zpro;

%% PLOT SNAPSHOT

ax = subplot(2,1,2); hold on
if (save_tag==0)
    if smooth_op == 1; pcolor(xmat./proZ,zmat./proZ,cmat);shading interp;
    elseif smooth_op == 2; pcolor(xsamp./proZ,zsamp./proZ,csamp);shading interp;
    else, patch(xrsh./proZ,zrsh./proZ,crsh,'EdgeColor','none'); 
    end; caxis([mV MV]); 
end

colormap(cmap);
t=colorbar('peer',gca);
set(get(t,'xlabel'),'String', 'D','Fontweight','bold')
annotation_plot('',touput,optplot,writcol,as,ax,2,2);

%% PLOT PROFILES

ax=subplot(2,1,1); hold on
for j=1:np
    %value on the profile +std
    plot([xp(j)*ones(Ny2/2,1)+Dprofil(1:Ny2/2,j)+Dstd(1:Ny2/2,j)-D01;xp(j)*ones(Ny2/2,1)+Dprofil(Ny2/2+1:end,j)+Dstd(Ny2/2+1:end,j)-D01],...
        Zprofil(:,j),'color',[0.7 0.7 0.7],'Linewidth',1)
    %value on the profile -std
    plot([xp(j)*ones(Ny2/2,1)+Dprofil(1:Ny2/2,j)-Dstd(1:Ny2/2,j)-D01;xp(j)*ones(Ny2/2,1)+Dprofil(Ny2/2+1:end,j)-Dstd(Ny2/2+1:end,j)-D01],...
        Zprofil(:,j),'color',[0.7 0.7 0.7],'Linewidth',1)
    %initial value
    plot([xp(j)*ones(Ny2/2,1)+D0profil(1:Ny2/2,j)-D01;xp(j)*ones(Ny2/2,1)+D0profil(Ny2/2+1:end,j)-D01],Zprofil(:,j),'k-')
    %value on the profile
    plot([xp(j)*ones(Ny2/2,1)+Dprofil(1:Ny2/2,j)-D01;xp(j)*ones(Ny2/2,1)+Dprofil(Ny2/2+1:end,j)-D01],Zprofil(:,j),'r-','Linewidth',1.5)
    %plot fault
    plot(fault.x./proZ,fault.z./proZ,'--','color',[0 0 0.7],'Linewidth',0.5)
end
xlim(xb); ylim(zb)
daspect([1 as 1])
colorbar
set(ax,'Xticklabel',[])
ylabel('Distance/R_0','FontSize', fs)
ax = gca;
ax.FontSize = fs;


%% SAVE FIGURES

%for illustrator
if (save_tag>=1) 
    
	%Save in eps format
    if (save_tag==1),export_fig(namefile,v_handle,'-eps'),end

    %Plot the screenshot
    ax=subplot(2,1,2); hold on
    if smooth_op == 1; pcolor(xmat./proZ,zmat./proZ,cmat);shading interp;
    elseif smooth_op == 2; pcolor(xsamp./proZ,zsamp./proZ,csamp);shading interp;
    else, patch(xrsh./proZ,zrsh./proZ,crsh,'EdgeColor','none'); 
    end; caxis([mV MV]); colormap(cmap);
    plot(fault.x./proZ,fault.z./proZ,'-w','Linewidth',0.5)
    t=colorbar('peer',gca);
    set(get(t,'xlabel'),'String', 'D','Fontweight','bold')
    xlim(xb); ylim(zb)
    daspect([1 as 1])

    %Comments
    annotation_plot('',touput,optplot,writcol,as,ax,2,2);
	
    %Save in matlab format
    saveas(v_handle,[namefile,'.fig'],'fig')

    %Save in png format
%     if (save_tag==1), export_fig(namefile,gcf,'-png','-r340');    end
    if (save_tag==1)
        ftmp=figure('Position',[200 500 1400 1200],'PaperOrientation','landscape');  
        if smooth_op == 1; pcolor(xmat./proZ,zmat./proZ,cmat);shading interp;...
                namefiletmp=[namefile,'_smooth1'];                
        elseif smooth_op == 2; pcolor(xsamp./proZ,zsamp./proZ,csamp);shading interp;...
                namefiletmp=[namefile,'_smooth2'];
        else, patch(xrsh./proZ,zrsh./proZ,crsh,'EdgeColor','none');...
                namefiletmp=namefile;               
        end; colormap(cmap); caxis([mV MV]);
        xlim(xb); ylim(zb);daspect([1 as 1])
        colorbar( 'off' )
        set(gca,'XTick',[], 'YTick', [])        
        export_fig(namefiletmp,gcf,'-png','-r340');
        close(figure(numfig+1));
    end

end

%% 
% =========================================================================
%                           LOG(D)-DISTANCE PLOT
% =========================================================================
%
%%
dy=fs/200;
xleg=zb(2)*0.64;
xleg2=zb(2)*0.68;

%colortable
ctab=[0 0 0;
    0.5 0 0.5;
%     0.8 0 0.8;
%     0.8 0.5 0.8;
    0.8 0.4 0.8;
%     1 0.7 0.7 ;
%     0.9 0.5 0.5;
    1 0.6 0.6 ;
    0.6 0 0;
    1 0 0;
    1 0.6 0;
    1 0.9 0;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
%     0.9 0.9 0.9;
%     0.8 0.8 0.8;
%     0.7 0.7 0.7;
    0.8 1 0;
%     0.6 0.8 0;
%     0.2 0.6 0;
    0.4 0.7 0;
    0 0.3 0.3;
    0 0 0.6;
    0 0.3 1;
    0 0.8 0.8;
%     0 1 1;
    0.8 1 1
	0 0 0];
% %launch figure
% figure; hold on
% plot(fault.x./proZ,fault.z./proZ,'-','color',[0 0 0.7],'Linewidth',0.5)
% plot(xp,zeros(size(xp)),'r*')

% find the point on the fault corresponding to the profiles
xfaultp=zeros(np,1);
zfaultp=zeros(np,1);
for i=1:np
    tmp=abs(fault.x./proZ-xp(i));
    idxf=find(tmp == min(tmp));
    xfaultp(i)=fault.x(idxf)./proZ;
    zfaultp(i)=fault.z(idxf)./proZ;
end
% plot(xfaultp,zfaultp,'go')

%Launch figure
v_handle =figure('Position',[200 500 900 900],'PaperOrientation','landscape');hold on
h=v_handle;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);

%title
subplot(2,2,1); hold on
title('upper left quadrant','FontSize', fs)
subplot(2,2,3); hold on
title('bottom left quadrant','FontSize', fs)
subplot(2,2,2); hold on
title('upper right quadrant','FontSize', fs)
subplot(2,2,4); hold on
title('bottom right quadrant','FontSize', fs)

%log log plot
i=0;k=1;

for j=1:np%2:np-1
    
    zp=Zprofil(:,j);d0p=D0pro(:,j);dp=Dpro(:,j);dsp=Dstdpro(:,j);

    if xp(j)<=0%-1*edge %xp(j)<=-2.5
    i=i+1;
    %upper left quadrant
    subplot(2,2,1); hold on
    idz=find(zp>=zfaultp(j));
    plot(zp(idz)-zfaultp(j),log10(dp(idz)),'o-','MarkerSize',4,'color',ctab(i,:),'Markerfacecolor',ctab(i,:));     %simulation value
    plot(zp(idz)-zfaultp(j),log10(d0p(idz)),'k-');         %initial value
    xlim([0 zb(2)]); 
    ylim([log10(0.095), log10(1)])
    plot(xleg,-dy*(j),'o','MarkerSize',6,'color',ctab(i,:),'Markerfacecolor',ctab(i,:))
    text(xleg2,-dy*(j),[' x = ', num2str(round(xp(j)*10)/10),' R_0'],'FontSize', fs,'Color','k');
%     set(gca,'XScale','Log');
    
    %bottom left quadrant 
    subplot(2,2,3); hold on
    idz=find(zp<zfaultp(j));
    plot(-zp(idz)+zfaultp(j),log10(dp(idz)),'o-','MarkerSize',4,'color',ctab(i,:),'Markerfacecolor',ctab(i,:));     % %simulation value
    plot(-zp(idz)+zfaultp(j),log10(d0p(idz)),'k-');         %initial value
    xlim([0 zb(2)]); 
    ylim([log10(0.095), log10(1)])
    plot(xleg,-dy*(j),'o','MarkerSize',6,'color',ctab(i,:),'Markerfacecolor',ctab(i,:))
    text(xleg2,-dy*(j),[' x = ', num2str(round(xp(j)*10)/10),' R_0'],'FontSize', fs,'Color','k');
%     set(gca,'XScale','Log');

    elseif xp(j)>0%=edge %xp(j)>=2.5
    i=i+1; k=k+1;
    
    %upper right quadrant
    subplot(2,2,2); hold on
    idz=find(zp>=zfaultp(j));
    plot(zp(idz)-zfaultp(j),log10(dp(idz)),'o-','MarkerSize',4,'color',ctab(i,:),'Markerfacecolor',ctab(i,:));     %simulation value
    plot(zp(idz)-zfaultp(j),log10(d0p(idz)),'k-');         %initial value
    xlim([0 zb(2)]); 
    ylim([log10(0.095), log10(1)])
    plot(xleg,-dy*(k-1),'o','MarkerSize',6,'color',ctab(i,:),'Markerfacecolor',ctab(i,:))
    text(xleg2,-dy*(k-1),[' x = ', num2str(round(xp(j)*10)/10),' R_0'],'FontSize', fs,'Color','k');
%     set(gca,'XScale','Log');
   
    %bottom right quadrant 
    subplot(2,2,4); hold on
    idz=find(zp<zfaultp(j));
    plot(-zp(idz)+zfaultp(j),log10(dp(idz)),'o-','MarkerSize',4,'color',ctab(i,:),'Markerfacecolor',ctab(i,:));     %%simulation value
    plot(-zp(idz)+zfaultp(j),log10(d0p(idz)),'k-');         %initial value
    xlim([0 zb(2)]); 
    ylim([log10(0.095), log10(1)])
    plot(xleg,-dy*(k-1),'o','MarkerSize',6,'color',ctab(i,:),'Markerfacecolor',ctab(i,:))
    text(xleg2,-dy*(k-1),[' x = ', num2str(round(xp(j)*10)/10),' R_0'],'FontSize', fs,'Color','k');
%     set(gca,'XScale','Log');
    
    end
end

if (save_tag>=1)
    export_fig([namef,'logDamage'],v_handle,'-eps'); 
    saveas(v_handle,[namef,'logDamage','.fig'],'fig');
end

