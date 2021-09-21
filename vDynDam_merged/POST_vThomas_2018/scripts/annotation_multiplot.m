% Script to annotate multiplot

% Marion Thomas last modified september 2018

%CALLS: 

%==========================================================================
%% 

colormap(ax,cmap);

%title
ax.FontSize = fs;
text((xb(2)-(xb(2)-xb(1))/1.8),(zb(2)-(zb(2)-zb(1))*0.1),...
    ['t = ', num2str(round(10*T(oo(o(i))+1))/10),'s'],'FontSize', fs-2,'Color',writcol);
    
% axis
xlim(xb); ylim(zb)
daspect([1 as 1])
ylabel('z/R_0','FontSize', fs)
set(ax,'YAxisLocation','right')
if (i==(j-1)*l+l || i==numel(o)), xticklabels('auto');xlabel('Distance/R_0','FontSize', fs);...
else, set(ax,'Xticklabel',[]);end
