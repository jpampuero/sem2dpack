% ************************************************************************
% This function plot the recorded variables such as stress, strain,
% particules acceleration or velocity, etc... it also plots the slip
% rate along the fault at the same time step.

% Marion Thomas, last modified September 2018

% CALLS: sem2d_snapshot_read.m; export_fig.m; sem2d_snapshot_plot_norm.m

%% ************************************************************************
function [] = plot_snapshots_nvar(datadir,grid,oo,varlist,T,fault,save_tag,namef,opt,cbar,wcol)

% datadir=datadirD;
% grid=gridD;
% oo=ooD;
% varlist=var;
% T=TD;
% fault=faultD;
% namef=namefoldD;
% opt=options;
% wcol='black';

%% INPUTS VARIABLES
%datadir:   directory of the data
%grid:      corresponding grid to plot the snapshot
%oo:        potential snapshots
%varlist        variables to be plotted
%T:         time
%fault:     simlations value on the fault
%save_tag:  saving option (1=save, 0=just plot, 2=save just .fig)
%namef:     folder
%cbar:      colormap for subplots
%writcol:	color for the writting

%options main plot
proZ = opt(16);         %process zone size
xb = opt(1:2)./proZ;    %min/max for x-axis
zb = opt(3:4)./proZ;    %min/max for y-axis
as = opt(5);            %aspect ratio between axis (x/y)
ym_v = opt(15);     	%max value for slip rate plot
tplot = opt(17);        %Time at which we want to plot the variable
smooth_op = opt(18);    %smoothing option for snapshot: 0 for NO, 1 for slight smoothing,
                        %2 for high smoothing like Kurama

%comments plot
optplot=[opt(1:4)./proZ,opt(6:14)];


%% PARAMETERS

%grid resolution
h=(fault.x(end)-fault.x(1))/(size(fault.x,1)-1);

%find snpashot to plot
o=find(abs(T-tplot)== min(abs(T-tplot)));
touput=num2str(round(10*T(oo(o)+1))/10);
texact=T(oo(o)+1);

%find number of variable to plot
n=numel(varlist);varnames='';
for j=1:n
    if j==n,varnames=[varnames,varlist{j}];else,...
            varnames=[varnames,varlist{j},' & '];end
end
disp(' ')
disp(['Snapshot at t=',num2str(tplot), 's for ', varnames])

%Launch figure
v_handle =figure('Position',[200 500 1400 1200],'PaperOrientation','landscape');  
f=v_handle;
objf=findobj('type','figure');numfig=length(objf);
set(f,'PaperOrientation','landscape');
set(f,'PaperUnits','normalized');
set(f,'PaperPosition', [0 0 1 1]);

%name files and folder
namefilemain='';
for j=1:n
    if j==n, namefilemain=[namefilemain,varlist{j}];else...
            namefilemain=[namefilemain,varlist{j},'_'];end
end
namef=[namef,'/',namefilemain,'/'];
if exist(namef) == 7; else; mkdir(namef); end
namefilemain=[namef,namefilemain,'_','snap',num2str(oo(o))];

%plot parameters
cmapstart=[1;cumsum(cbar(end-n+1:end-1,1))+1];
cmapend=cumsum(cbar(end-n+1:end,1));
cmapval=[cmapstart cmapend];
mV=zeros(1,n);
MV=zeros(1,n);

%% DOWNLOAD

for j=1:n

    %download variable
    i=1;var=varlist{j};
    var2plot

    %min/max
    if size(opt,2) == 18+2*n
        mV(j) = opt(18+2*(j-1)+1);          %maximum value for the 1st variable
        MV(j) = opt(18+2*(j-1)+2);          %minimum value for the 1st variable
    else
        MV(j)=max(max(max(abs(v))));
        mV(j)=-MV(j);%min(min(min(v)));
    end
     
    %% PLOT VARIABLES

    %Subfigures
    ax=subplot(n,1,j); hold on

    if (save_tag==0)
        
       %Matrix rearrangement
       if smooth_op >= 1; data_matrix_sorting; end
       
        %plotting
        if smooth_op == 1; pcolor(xmat./proZ,zmat./proZ,cmat);shading interp;
        elseif smooth_op == 2; pcolor(xsamp./proZ,zsamp./proZ,csamp);shading interp;
%        % draw contourlines
%         contourresolution_l = linspace(mV,MV,20);
%         [~, hpl] = contourf(xsamp./proZ,zsamp./proZ,csamp,contourresolution_l,...
%             '-', 'LineWidth',0.5,'LineColor',[1 1 1]*0.35);
%          set(hpl,'Fill','off');
        else, sem2d_snapshot_plot_norm(v,grid,proZ,[mV MV]);
        end
    end
    
    %colorbar/annotation
    colormap(ax,cbar(cmapval(j,1):cmapval(j,2),:)); caxis([mV(j) MV(j)]);
    t=colorbar('peer',gca); 
    set(get(t,'xlabel'),'String', barvar,'Fontweight','bold')
    if j==1; annotation_plot(vname,touput,optplot,writcol,as,ax,n,j); ...
    else, annotation_plot2(vname,optplot,writcol,as,ax,n,j); end

    %fault
    plot(fault.x./proZ,fault.z./proZ,'w','Linewidth',0.5)

    %normalized slip rate
    id=find(abs(fault.t-texact)== min(abs(fault.t-texact)));
%     plot(fault.x./proZ,fault.v(:,id)./ym_v,'Color',colorL,'Linewidth',1)
    plot(fault.x./proZ,fault.v(:,id)./ym_v*zb(2),'Color',colorL,'Linewidth',1)

end

%% SAVE FIGURES for illustrator

if (save_tag>=1) 
    
    %Save in eps format
    if (save_tag==1),export_fig(namefilemain,v_handle,'-eps'),end

    for j=1:n
    
    %download variable
    i=1;var=varlist{j};
    var2plot

	%Matrix rearrangement
	if smooth_op >= 1; data_matrix_sorting; end
        
    % Plot the variable
    subplot(n,1,j)
    if smooth_op == 1; pcolor(xmat./proZ,zmat./proZ,cmat);shading interp;
    elseif smooth_op == 2; pcolor(xsamp./proZ,zsamp./proZ,csamp);shading interp;
    else, sem2d_snapshot_plot_norm(v,grid,proZ,[mV MV]);
    end
    colormap(ax,cbar(cmapval(j,1):cmapval(j,2),:)); caxis([mV(j) MV(j)]);
  
    %Comments
    t=colorbar('peer',gca); 
    set(get(t,'xlabel'),'String', barvar,'Fontweight','bold')
    if j==1; annotation_plot(vname,touput,optplot,writcol,as,ax,n,j); ...
    else, annotation_plot2(vname,optplot,writcol,as,ax,n,j); end

    %fault
    plot(fault.x./proZ,fault.z./proZ,'w','Linewidth',0.5)

    %normalized slip rate
    plot(fault.x./proZ,fault.v(:,id)./ym_v,'Color',colorL,'Linewidth',1)

    %Save in png format
    if (save_tag==1)
        ftmp=figure('Position',[200 500 1400 1200],'PaperOrientation','landscape');  
        if smooth_op == 1; pcolor(xmat./proZ,zmat./proZ,cmat);shading interp;...
                namefiletmp=[namef,varlist{j},'_','snap',num2str(oo(o)),'_smooth1'];
        elseif smooth_op == 2; pcolor(xsamp./proZ,zsamp./proZ,csamp);shading interp;...
                namefiletmp=[namef,varlist{j},'_','snap',num2str(oo(o)),'_smooth2'];
        else, sem2d_snapshot_plot_norm(v,grid,proZ,[mV MV]);...
                namefiletmp=[namef,varlist{j},'_','snap',num2str(oo(o))];
        end
        colormap(gca,cbar(cmapval(j,1):cmapval(j,2),:)); caxis([mV(j) MV(j)]);
        xlim(xb); ylim(zb);daspect([1 as 1])
        colorbar( 'off' )
        set(gca,'XTick',[], 'YTick', [])        
        export_fig(namefiletmp,gcf,'-png','-r340')
        close(figure(numfig+1));
    end
    end
    
    %Save in matlab format
    saveas(v_handle,[namefilemain,'.fig'],'fig')
    
end
