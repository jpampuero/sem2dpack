% ************************************************************************
% This function plot useful information on the plot

% Marion Thomas, last modified May 2017

% CALLS: sem2d_snapshot_read.m; export_fig.m; sem2d_snapshot_plot_km.m

%% ************************************************************************
function [] = annotation_plot2(vn,opt,writcol,as,ax,nv,i)

%% INPUTS VARIABLES

%vn:        Name of the plotted variable
%touput:    time of the snapshot
%writcol:	color for the writting

%options main plot
xb = opt(1:2);  %min/max for x-axis
yb = opt(3:4);  %min/max for y-axis
fs = opt(5); % font size

%% 
dy=fs/180;
dx=fs/100;

%Comments
text(xb(2)-((xb(2)-xb(1))*1.5*dx),yb(1)+((yb(2)-yb(1))*dy),['$$' vn '$$'],'FontSize', fs,'interpreter','latex','Color',writcol);

%Arrows
text(0,0.125,'$$\rightarrow$$','interpreter','latex','color','w','fontsize',21)
text(0,-0.125,'$$\leftarrow$$','interpreter','latex','color','w','fontsize',21)

% axis
xlim(xb); ylim(yb)
daspect([1 as 1])
ylabel('Distance/R_0','FontSize', fs)
% set(ax,'YAxisLocation','right')
if (i==nv), xticklabels('auto');xlabel('Distance/R_0','FontSize', fs);...
else, set(ax,'Xticklabel',[]);end

ax = gca;
ax.FontSize = fs;

