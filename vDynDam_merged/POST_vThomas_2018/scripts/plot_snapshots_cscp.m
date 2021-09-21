% *************************************************************************
% This function reads plot the change in Cs and Cp velocity due to damage

% Marion Thomas, last modified September 2018

% CALLS: sem2d_snapshot_read.m; export_fig.m; sem2d_snapshot_plot_km.m

%% ************************************************************************
function [] = plot_snapshots_cscp(datadir,grid,oo,T,fault,save_tag,namef,opt,writcol,cmap)

%% INPUTS VARIABLES
%datadir:   directory of the data
%grid:      corresponding grid to plot the snapshot
%oo:        potential snapshots
%var1,var2  variables to be plotted
%T:         time
%fault: 	simlations value on the fault
%save_tag:  saving option (1=save, 0=just plot, 2=save just .fig)
%namef:     folder
%cmap:      colormap for subplots
%writcol:	color for the writting

%options main plot
proZ = opt(15);         %process zone size
xb = opt(1:2)./proZ;    %min/max for x-axis
zb = opt(3:4)./proZ;    %min/max for y-axis
as = opt(5);            %aspect ratio between axis (x/y)
tplot = opt(16);        %Time at which we want to plot the variable
DThres = opt(17);       %Threshold after which the contribution of interation between cracks takes over
smooth_op = opt(18);    %smoothing option for snapshot: 1 for yes 0 for no.

%comments plot
optplot=[opt(1:4)./proZ,opt(6:14)];
disp(' ')
disp('Reduction in S- and P-waves are plotted')

%% PARAMETERS

%grid resolution
h=(fault.x(end)-fault.x(1))/(size(fault.x,1)-1);

%find snpashot to plot
o=find(abs(T-tplot)== min(abs(T-tplot)));
touput=num2str(round(10*T(oo(o)+1))/10);

%Launch figure
v_handle =figure('Position',[200 500 1600 750],'PaperOrientation','landscape');  
f=v_handle;
objf=findobj('type','figure');numfig=length(objf);
set(f,'PaperOrientation','landscape');
set(f,'PaperUnits','normalized');
set(f,'PaperPosition', [0 0 1 1]);

%name files and folder
namef=[namef,'/','cs_cp_variation/'];
if exist(namef) == 7; else; mkdir(namef); end
namefile=[namef,'cs_cp_variation_snap',num2str(oo(o))];

%% DOWNLOAD

%downlaod
dyD = sem2d_snapshot_read('dyD',oo(o),datadir);

%Find the variable to plot
v1=dyD.cs;
cs0=max(max(max(v1)));
v1=100*(cs0-v1)/cs0;%max(max(max(v1))); 
vn1='\%\ of\ reduction\ in\ c_s';%'\%\ of\ reduction\ in\ S-waves\ speed';
v2=dyD.cp;
cp0=max(max(max(v2)));
v2=100*(cp0-v2)/cp0;%max(max(max(v2)));
% vn2='\%\ of\ reductio\n in\ P-waves\ speed';
vn2='\%\ of\ reduction\ in\ c_p';%'\%\ of\ reduction\ in\ S-waves\ speed';

% min/max
if size(opt,2) == 22
    mV1 = opt(19);          %maximum value for the main variable
    MV1 = opt(20);          %minimum value for the main variable
    mV2 = opt(21);          %minimum value for the main variable
    MV2 = opt(22);          %maximum value for the main variable
else
    MV1=ceil(max(max(max(v1))));
    mV1=floor(min(min(min(v1))));
    MV2=ceil(max(max(max(v2))));
    mV2=floor(min(min(min(v2))));
end

%% CORRECTION FOR CP CS

D=dyD.D;
idD=find(D<DThres);
v1(idD)= 0;
v2(idD)= 0;

%% PLOT Cs

%subplot
ax=subplot(2,1,1); hold on

%Matrix rearrangement
if smooth_op >= 1; v=v1; data_matrix_sorting; end

%plotting
if (save_tag==0)
    if smooth_op == 1; pcolor(xmat./proZ,zmat./proZ,cmat);shading interp;
    elseif smooth_op == 2; pcolor(xsamp./proZ,zsamp./proZ,csamp);shading interp;
    else, sem2d_snapshot_plot_norm(v1,grid,proZ,[mV1 MV1]);
    end
end
    
%colorbar/annotation
colormap(cmap); caxis([mV1 MV1]);
t=colorbar('peer',gca);
set(get(t,'xlabel'),'String', '%','Fontweight','bold')
annotation_plot(vn1,touput,optplot,writcol,as,ax,2,1);
xlim(xb); ylim(zb)
daspect([1 as 1])

%fault
plot(fault.x./proZ,fault.z./proZ,'-k','Linewidth',0.5)

%% PLOT Cp

%subplot
ax=subplot(2,1,2); hold on

%Matrix rearrangement
if smooth_op >= 1; v=v2; data_matrix_sorting; end

%plotting
if (save_tag==0)
    if smooth_op == 1; pcolor(xmat./proZ,zmat./proZ,cmat);shading interp;
    elseif smooth_op == 2; pcolor(xsamp./proZ,zsamp./proZ,csamp);shading interp;
    else, sem2d_snapshot_plot_norm(v2,grid,proZ,[mV2 MV2]);
    end
end
    
%colorbar/annotation
colormap(cmap); caxis([mV2 MV2]);
t=colorbar('peer',gca);
set(get(t,'xlabel'),'String', '%','Fontweight','bold')
annotation_plot2(vn2,optplot,writcol,as,ax,2,2);
xlim(xb); ylim(zb)
daspect([1 as 1])

%fault
plot(fault.x./proZ,fault.z./proZ,'-k','Linewidth',0.5)

%% SAVE FIGURES for illustrator

if (save_tag>=1) 
    
    %Save in eps format
    if (save_tag==1),export_fig(namefile,v_handle,'-eps'),end

    % PLOT Cs
    subplot(2,1,1);hold on
    if smooth_op >= 1; v=v1; data_matrix_sorting; end
    if smooth_op == 1; pcolor(xmat./proZ,zmat./proZ,cmat);shading interp;
    elseif smooth_op == 2; pcolor(xsamp./proZ,zsamp./proZ,csamp);shading interp;
    else, sem2d_snapshot_plot_norm(v1,grid,proZ,[mV1 MV1]);
    end
	plot(fault.x./proZ,fault.z./proZ,'-k','Linewidth',0.5)
    annotation_plot(vn1,touput,optplot,writcol,as,ax,2,1);
    
	%Save in png format
	if (save_tag==1)
        ftmp=figure('Position',[200 500 1400 1200],'PaperOrientation','landscape');  
        if smooth_op == 1; pcolor(xmat./proZ,zmat./proZ,cmat);shading interp;...
                namefiletmp=[namef,'cs_snap',num2str(oo(o)),'_smooth1'];
        elseif smooth_op == 2; pcolor(xsamp./proZ,zsamp./proZ,csamp);shading interp;...
                namefiletmp=[namef,'cs_snap',num2str(oo(o)),'_smooth2'];
        else, sem2d_snapshot_plot_norm(v1,grid,proZ,[mV1 MV1]);...
                namefiletmp=[namef,'cs_snap',num2str(oo(o))];
        end
        colormap(cmap); caxis([mV1 MV1]);colorbar( 'off' )
        xlim(xb); ylim(zb);daspect([1 as 1])
        set(gca,'XTick',[], 'YTick', [])                
        export_fig(namefiletmp,gcf,'-png','-r340')
        close(figure(numfig+1));
	end

    % PLOT Cp
    subplot(2,1,2); hold on
    if smooth_op >= 1; v=v1; data_matrix_sorting; end
    if smooth_op == 1; pcolor(xmat./proZ,zmat./proZ,cmat);shading interp;
    elseif smooth_op == 2; pcolor(xsamp./proZ,zsamp./proZ,csamp);shading interp;
    else, sem2d_snapshot_plot_norm(v2,grid,proZ,[mV2 MV2]);
    end
	plot(fault.x./proZ,fault.z./proZ,'-k','Linewidth',0.5)
    annotation_plot2(vn2,optplot,writcol,as,ax,2,2);    

	%Save in png format
	if (save_tag==1)
        ftmp=figure('Position',[200 500 1400 1200],'PaperOrientation','landscape');  
        if smooth_op == 1; pcolor(xmat./proZ,zmat./proZ,cmat);shading interp;...
                namefiletmp=[namef,'cp_snap',num2str(oo(o)),'_smooth1'];
        elseif smooth_op == 2; pcolor(xsamp./proZ,zsamp./proZ,csamp);shading interp;...
                namefiletmp=[namef,'cp_snap',num2str(oo(o)),'_smooth2'];
        else, sem2d_snapshot_plot_norm(v2,grid,proZ,[mV2 MV2]);...
                namefiletmp=[namef,'cp_snap',num2str(oo(o))];
        end
        colormap(cmap); caxis([mV1 MV1]);colorbar( 'off' )
        xlim(xb); ylim(zb);daspect([1 as 1])
        set(gca,'XTick',[], 'YTick', [])                
        export_fig(namefiletmp,gcf,'-png','-r340')
        close(figure(numfig+1));
	end
        
        
    %Save in matlab format
    saveas(v_handle,[namefile,'.fig'],'fig')
    
end
