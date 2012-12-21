% Plot2dSnapshot(field,coord,indx,vsat)
%
% INPUT		coord	2D coordinates of GLL nodes (npoin,2)
%			field	field data to be plotted
%			indx	cell info, output from Init2dSnapshot
%			vsat	color scale [min max], default none
%
% OUTPUT	h	patch handle 
%
function h=Plot2dSnapshot(field,coord,indx,vsat)

set(gca,'DefaultPatchEdgeColor','none');
h=patch('Vertices',coord,'Faces',indx,'FaceVertexCData',field(:),'FaceColor','interp');
axis equal tight
title('SEM2D snapshot')
xlabel('X')
ylabel('Z')
if exist('vsat','var'), if ~isempty(vsat), set(gca,'CLim',vsat); end, end
colorbar('vert')
