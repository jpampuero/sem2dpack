% SEM2D_SNAPSHOT_PLOT color-plots a snapshot output from SEM2DPACK
%
% SYNTAX	sem2d_snapshot_plot(field,grid,vsat)
%
% INPUT		field	snapshot data, single field 
%			in element-wise or node-wise format
%		grid	spectral element grid (see SEM2D_READ_SPECGRID)
%		vsat	['minmax'] color scale
%			'minmax' = [min(field) max(field)]
%			'maxabs' = [-1 1]*max(abs(field))
%			[lo hi]  = prescribed numeric bounds
%
function sem2d_snapshot_plot(field,grid,vsat)

set(gca,'DefaultPatchEdgeColor','none');

if ndims(field) == 3
  [indx,coord_new] = initialize(grid.ibool,grid.coord);
  h=patch('Vertices',coord_new,'Faces',indx,'FaceVertexCData',field(:),'FaceColor','interp');
else
  indx = initialize(grid.ibool);
  h=patch('Vertices',grid.coord,'Faces',indx,'FaceVertexCData',field(:),'FaceColor','interp');
end

% to plot also the grid:
%hold on; plot(grid.coord(:,1),grid.coord(:,2),'k+'); hold off

axis equal tight
box on
title('SEM2D snapshot')
xlabel('X')
ylabel('Z')

if ~exist('vsat','var') | isempty(vsat), vsat = 'minmax'; end
if isstr(vsat)
  switch lower(vsat)
    case 'maxabs',
      vsat = [-1 1] * max(abs(field(:)));
    case 'minmax',
      vsat = [ min(field(:)) max(field(:)) ];
    otherwise,
      error('Unknown value for argument vsat')
  end
end
vsat = vsat(1:2);
vsat = vsat(:)';
if vsat(1)==vsat(2), vsat = vsat+[-1 1]; end
set(gca,'CLim',vsat)
colorbar('vert')


%------------------------------------------------------------------
% indx = initialize(iglob)
% [indx,coord_new] = initialize(iglob,coord)
%
% PURPOSE	Initializes grid data for Plot2dSnapshot
%
% INPUT		iglob(ngll,ngll,nel) local to global numbering map
%		coord(npoin,2) 	coordinates of the nodes
%
% OUTPUT	indx(:,4) 	vertices of each GLL cell
%		coord_new(:,2) 	reshaped coord 
%
% NOTE		coord and coord_new are required only to plot element-wise fields
%		(stress, strain, material properties)
%
function [indx,coord_new] = initialize(iglob,coord)

[NGLL,NGLL,NEL] = size(iglob);
indx = zeros(NEL*(NGLL-1)^2, 4);

% node-wise storage
if nargin<2

  ip=0;
  for e=1:NEL,
    for i=1:NGLL-1,
    for j=1:NGLL-1,
      ip = ip+1;
      indx(ip,:) = [iglob(i,j,e) iglob(i+1,j,e) iglob(i+1,j+1,e) iglob(i,j+1,e)];
    end
    end
  end

% element-wise storage
else
 
  if nargout<2, error('Too few outputs'), end

  coord_new = coord(iglob(:),1:2);

  iglob0 = reshape( (1:NGLL*NGLL)', NGLL,NGLL);
  indx0 = zeros((NGLL-1)^2, 4);
  ip=0;
  for i=1:NGLL-1,
  for j=1:NGLL-1,
    ip = ip+1;
    indx0(ip,:) = [iglob0(i,j) iglob0(i+1,j) iglob0(i+1,j+1) iglob0(i,j+1)];
  end
  end

  blocksize = (NGLL-1)^2;
  for e=1:NEL,
    ip = (e-1)*blocksize;
    indx(ip+1:ip+blocksize,:) = indx0 + (e-1)*NGLL^2;
  end

end
