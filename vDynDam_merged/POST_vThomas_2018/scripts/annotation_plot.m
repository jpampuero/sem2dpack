% ************************************************************************
% This function plot useful information on the plot

% Marion Thomas, last modified May 2017

% CALLS: sem2d_snapshot_read.m; export_fig.m; sem2d_snapshot_plot_km.m

%% ************************************************************************
function [] = annotation_plot(vn,touput,opt,writcol,as,ax,nv,i)

%% INPUTS VARIABLES

%vn:        Name of the plotted variable
%touput:    time of the snapshot
%writcol:	color for the writting

%options main plot
xb = opt(1:2);  %min/max for x-axis
yb = opt(3:4);  %min/max for y-axis
fs = opt(5); % font size

D01 = round(10*opt(6))/10;    %initial damage for lower part
cs1 = round(10*opt(7)/1e3)/10;   %initial cs for lower part
cp1 = round(10*opt(8)/1e3)/10;   %initial cp for lower part
rho1 = round(10*opt(9))/10;	%initial rho for lower part
D02 = round(10*opt(10))/10;   %initial damage for upper part
cs2 = round(10*opt(11)/1e3)/10;   %initial cs for upper part
cp2 = round(10*opt(12)/1e3)/10;   %initial cp for upper part
rho2 = round(10*opt(13))/10;    %initial rho for upper part

%% 
% dy=fs/230;
% dx=fs/140;
dy=fs/180;
dx=fs/100;

%Comments
text((xb(2)-(xb(2)-xb(1))/2)-3*dx,yb(2)-((yb(2)-yb(1))*dy),[' t = ', touput,' s'],'FontSize', fs,'interpreter','latex','Color',writcol);
text(xb(2)-((xb(2)-xb(1))*1.5*dx),yb(1)+((yb(2)-yb(1))*dy),['$$' vn '$$'],'FontSize', fs,'interpreter','latex','Color',writcol);

%properties lower part
text(xb(1)+((xb(2)-xb(1))*0.01),yb(1)+((yb(2)-yb(1))*dy),['$$D_0 = $$', num2str(D01)],'FontSize', fs,'interpreter','latex','Color',writcol);
text(xb(1)+((xb(2)-xb(1))*0.01),yb(1)+((yb(2)-yb(1))*dy*2),['$$c_s = $$', num2str(cs1),' km/s'],'FontSize', fs,'interpreter','latex','Color',writcol);
text(xb(1)+((xb(2)-xb(1))*0.01),yb(1)+((yb(2)-yb(1))*dy*3),['$$c_p = $$', num2str(cp1),' km/s'],'FontSize', fs,'interpreter','latex','Color',writcol);
text(xb(1)+((xb(2)-xb(1))*0.01),yb(1)+((yb(2)-yb(1))*dy*4),['$$\rho = $$', num2str(rho1),' kg/m','$$^3$$'],'FontSize', fs,'interpreter','latex','Color',writcol);

%properties upper part
coefy=1;
if rho1~=rho2
    dyy=dy*coefy;coefy=coefy+1;
    text(xb(1)+((xb(2)-xb(1))*0.01),yb(2)-((yb(2)-yb(1))*dy),['$$\rho = $$', num2str(rho2),' kg/m','$$^3$$'],'FontSize', fs,'interpreter','latex','Color',writcol);
end
if cp1~=cp2
    dyy=dy*coefy;coefy=coefy+1;
    text(xb(1)+((xb(2)-xb(1))*0.01),yb(2)-((yb(2)-yb(1))*dyy),['$$c_p = $$', num2str(cp2),' km/s'],'FontSize', fs,'interpreter','latex','Color',writcol);
end
if cs1~=cs2
    dyy=dy*coefy;coefy=coefy+1;
    text(xb(1)+((xb(2)-xb(1))*0.01),yb(2)-((yb(2)-yb(1))*dyy),['$$c_s = $$', num2str(cs2),' km/s'],'FontSize', fs,'interpreter','latex','Color',writcol);
end
if D01~=D02
    dyy=dy*coefy;coefy=coefy+1;
    text(xb(1)+((xb(2)-xb(1))*0.01),yb(2)-((yb(2)-yb(1))*dyy),['$$D_0 = $$', num2str(D02)],'FontSize', fs,'interpreter','latex','Color',writcol);
end

%Arrows
text((xb(2)-(xb(2)-xb(1))/2)-1*dx,0.125,'$$\rightarrow$$','interpreter','latex','color','w','fontsize',21)
text((xb(2)-(xb(2)-xb(1))/2)-1*dx,-0.125,'$$\leftarrow$$','interpreter','latex','color','w','fontsize',21)

% axis
xlim(xb); ylim(yb)
daspect([1 as 1])
ylabel('Distance/R_0','FontSize', fs)
% set(ax,'YAxisLocation','right')
if (i==nv), xticklabels('auto');xlabel('Distance/R_0','FontSize', fs);...
else, set(ax,'Xticklabel',[]);end

ax = gca;
ax.FontSize = fs;

