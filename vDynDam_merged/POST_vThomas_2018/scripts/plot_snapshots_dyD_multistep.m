% ************************************************************************
% This function plot n snapshosts of recorded dynamic damage variables such
% as damage, change in cs or cp, dl/dt, KI, etc... 

% Marion Thomas, last modified September 2018

% CALLS: sem2d_snapshot_read.m; export_fig.m; sem2d_snapshot_plot_km.m

%% ************************************************************************
function [] = plot_snapshots_dyD_multistep(datadir,grid,n,oo,var,T,fault,save_tag,namef,opt,cmap,wcol)


%% INPUTS VARIABLES
%datadir:   directory of the data
%grid:      corresponding grid to plot the snapshot
%n:         number of snapshots
%oo:        potential snapshots
%var:       variable to be plotted
%T:         time for the snapshots
%save_tag:  saving option (1=save, 0=just plot, 2=save just .fig)
%namef      Name of the subfolder for output
%cmap:      colormap for subplots
%writcol:	color for the writting

%options 
proZ = opt(8);         %process zone size
xb = opt(1:2)./proZ;    %min/max for x-axis
zb = opt(3:4)./proZ;    %min/max for y-axis
as = opt(5);            %aspect ratio between axis (x/y)
ym_v = opt(7);          %max value for slip rate plot
tstart = opt(9); 	    %Time at which we want to start plotting the variable
tend = opt(10);        %Time at which we want to end plotting the variable
smooth_op = opt(11);    %smoothing option for snapshot: 0 for NO, 1 for slight smoothing,
                        %2 for high smoothing like Kurama

%comments plot
fs=opt(6);

%% PLOT PARAMETERS

%grid resolution
h=(fault.x(end)-fault.x(1))/(size(fault.x,1)-1);

%find snpashot to plot
o_end=find(abs(T-tend)== min(abs(T-tend)));
o_sta=find(abs(T-tstart)== min(abs(T-tstart)));

%find snpashot to plot
step=round((o_end-o_sta)/(n-1));
o=sort(o_end:-step:o_sta);
n=numel(o)+1;

%number of subplot
if n<=2
    l=n; c=1; mat_idx=1:n;
else
    l = ceil(n/2); ceil(sqrt(n)); 
    c = ceil(n/l);
    %sort index for subplot
    mat_idx=reshape([1:numel(o),numel(o)+1:c*l],[c,l])';
end

%Launch figure
v_handle =figure('Position',[200 500 1600 1000],'PaperOrientation','landscape');  
f=v_handle;
objf=findobj('type','figure');numfig=length(objf);
set(f,'PaperOrientation','landscape');
set(f,'PaperUnits','normalized');
set(f,'PaperPosition', [0 0 1 1]);

%name files and folder
namef=[namef,'/',var,'_multi_snapshots/'];
if exist(namef) == 7; else; mkdir(namef); end
namefilemain=[namef,var,'_multi_snapshots'];

%comment
disp(' ')
disp([num2str(n-1),' snapshots of ', var, ' from t=',num2str(tstart), 's to ', num2str(tend), 's'])

%plot parameters
mV=zeros(1,n-1);
MV=zeros(1,n-1);

%% DOWNLOAD

for j=1:c    
    for k=1:l        
    i=(j-1)*l+k;
    
    if i<=numel(o)

    %downlaod
    dyD = sem2d_snapshot_read('dyD',oo(o(i)),datadir);
    var2plot

    % min/max
    if mat_idx(i)==1
    if size(opt,2) >= 13
        mV = opt(12);          %maximum value for the variable
        MV = opt(13);          %minimum value for the variable
    else
        MV=max(max(max(v)));
        mV=min(min(min(v)));
    end
    end
    
%% PLOT VARIABLES

    %Subfigures
    ax = subplot(l,c,mat_idx(i)); hold on;
	if (save_tag==0)
     
	%Matrix rearrangement
        if smooth_op >= 1; data_matrix_sorting; end
       
	%plotting
        if smooth_op == 1; pcolor(xmat./proZ,zmat./proZ,cmat);shading interp;
        elseif smooth_op == 2; pcolor(xsamp./proZ,zsamp./proZ,csamp);shading interp;
        else, sem2d_snapshot_plot_norm(v,grid,proZ,[mV MV]);
        end
    end

	%fault
    plot(fault.x./proZ,fault.z./proZ,'w','Linewidth',0.5)

    %normalized slip rate
    temp=abs(fault.t-T(oo(o(i))+1));
    id=find(temp == min(temp));
    plot(fault.x./proZ,fault.v(:,id)./ym_v,'Color',colorL,'Linewidth',0.5)

    %comment plot
    annotation_multiplot
    end
    
    end
end

%plot for the colorbar
ax = subplot(max(l),max(c),max(l)*max(c)); hold on
xlim(xb); ylim(zb); daspect([1 as 1])
colormap(ax,cmap);caxis([mV MV]);
colorbar('southoutside');

%% SAVE FIGURES for illustrator

%for illustrator
if (save_tag>=1) 

    %Save in eps format
    if (save_tag==1),export_fig(namefilemain,v_handle,'-eps'),end

    for j=1:c
        for k=1:l
        i=(j-1)*l+k;

        if i<=numel(o)

        %downlaod variable
        dyD = sem2d_snapshot_read('dyD',oo(o(i)),datadir);
        var2plot

        %Matrix rearrangement
        if smooth_op >= 1; data_matrix_sorting; end

        %plot the variable
        ax = subplot(l,c,mat_idx(i)); hold on;
        if smooth_op == 1; pcolor(xmat./proZ,zmat./proZ,cmat);shading interp;
        elseif smooth_op == 2; pcolor(xsamp./proZ,zsamp./proZ,csamp);shading interp;
        else, sem2d_snapshot_plot_norm(v,grid,proZ,[mV MV]);
        end

        %fault
        plot(fault.x./proZ,fault.z./proZ,'w','Linewidth',0.5)

        %normalized slip rate
        temp=abs(fault.t-T(oo(o(i))+1));
        id=find(temp == min(temp));
        plot(fault.x./proZ,fault.v(:,id)./ym_v,'Color',colorL,'Linewidth',0.5)

        %comment plot
        annotation_multiplot
        
         %Save in png format
        if (save_tag==1)
            ftmp=figure('Position',[200 500 1400 1200],'PaperOrientation','landscape');  
            if smooth_op == 1; pcolor(xmat./proZ,zmat./proZ,cmat);shading interp;...
                    namefiletmp=[namef,var,'_','snap',num2str(oo(o(i))+1),'_smooth1'];
            elseif smooth_op == 2; pcolor(xsamp./proZ,zsamp./proZ,csamp);shading interp;...
                    namefiletmp=[namef,var,'_','snap',num2str(oo(o(i))+1),'_smooth2'];
            else, sem2d_snapshot_plot_norm(v,grid,proZ,[mV MV]);...
                    namefiletmp=[namef,var,'_','snap',num2str(oo(o(i))+1)];
            end
            colormap(gca,cmap); caxis([mV MV]);
            xlim(xb); ylim(zb);daspect([1 as 1])
            colorbar( 'off' )
            set(gca,'XTick',[], 'YTick', [])        
            export_fig(namefiletmp,gcf,'-png','-r340')
            close(figure(numfig+1));
        end
        
        end
        end
        
    end
    
	saveas(v_handle,[namefilemain,'.fig'],'fig')

    %     %Save in png format
%     if (save_tag==1), export_fig(namefilemain,gcf,'-png','-r340');end

%     for j=1:c
%         for k=1:l
%         i=(j-1)*l+k;
% 
%         if i<=numel(o)
% 
%         %normalized slip rate/fault
%         ax = subplot(l,c,mat_idx(i)); hold on
%         temp=abs(fault.t-T(oo(o(i))+1));
%         id=find(temp == min(temp));
%         plot(fault.x./proZ,fault.v(:,id)./ym_v,'r','Linewidth',1)
%         plot(fault.x./proZ,fault.z./proZ,'w','Linewidth',1)
% 
%         %comment plot
%         annotation_multiplot
%         end
% 
%         end
% 	end

end

end


%% Colorbar

% v_handle2 =figure('Position',[200 500 1600 750],'PaperOrientation','landscape');  
% f2=v_handle2; hold on
% set(f2,'PaperOrientation','landscape');
% set(f2,'PaperUnits','normalized');
% set(f2,'PaperPosition', [0 0 1 1]);
% namefile2=[namef,var,'_last_snapshots'];
% 
% %downlaod
% dyD = sem2d_snapshot_read('dyD',oo(o(end)),datadir);
% i=numel(o);var2plot
%    
% %plot
% if (save_tag==0),sem2d_snapshot_plot_norm(v,grid,proZ,[mV MV]);end
% plot(fault.x./proZ,fault.z./proZ,'w','Linewidth',1)
%     
% %normalized slip rate
% temp=abs(fault.t-T(oo(o(end))+1));
% id=find(temp == min(temp));
% plot(fault.x./proZ,fault.v(:,id)./ym_v,'r','Linewidth',1)
% plot(fault.x./proZ,fault.z./proZ,'w','Linewidth',0.5)
% 
% %comment plot
% colormap(cmap); colorbar('northoutside')
% xlim(xb); ylim(yb); daspect([1 as 1])
% ylabel('z/R_0','FontSize', fs); 
% set(ax,'YAxisLocation','right');xlabel('Distance/R_0','FontSize', fs);
% 
% if (save_tag>=1) 
%     
%     if (save_tag==1),export_fig(namefile2,v_handle2,'-eps');end
% 
%     % Plot the 3rd variable
%     sem2d_snapshot_plot_norm(v,grid,proZ,[mV MV]);    
%     xlim(xb); ylim(yb)
%     daspect([1 as 1])
%     colormap(cmap);
% 
%     if (save_tag==1), export_fig(namefile2,gcf,'-png','-r340');end
%     
%     %normalized slip rate / fault
%     plot(fault.x./proZ,fault.z./proZ,'w','Linewidth',1)
%     temp=abs(fault.t-T(oo(o(end))+1));
%     id=find(temp == min(temp));
%     plot(fault.x./proZ,fault.v(:,id)./ym_v,'r','Linewidth',1)
%     plot(fault.x./proZ,fault.z./proZ,'w','Linewidth',0.5)
%     
% 	saveas(v_handle2,[namefile2,'.fig'],'fig')
% 
% end

