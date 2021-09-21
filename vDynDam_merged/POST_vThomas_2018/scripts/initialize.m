%------------------------------------------------------------------
% indx = initialize(iglob)
% [indx,coord_new] = initialize(iglob,coord,target_box)
%
% PURPOSE	Initializes grid data for Plot2dSnapshot
%
% INPUT		iglob(ngll,ngll,nel) local to global numbering map
%		coord(npoin,2) 	coordinates of the nodes
%		target_box(4)	[xmin xmax ymin ymax]
%
% OUTPUT	indx(:,4) 	vertices of each GLL cell
%		coord_new(:,2) 	reshaped coord 
%
% NOTE		coord and coord_new are required only to plot element-wise fields
%		(stress, strain, material properties)
%
function [indx,coord_new] = initialize(iglob,coord,target_box)

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
  if exist('target_box','var') && ~isempty(target_box),
    kmid = floor(NGLL/2);
    xmin = target_box(1);
    xmax = target_box(2);
    ymin = target_box(3);
    ymax = target_box(4);
    ip = 0;
    for e=1:NEL,
      midpoint = coord( iglob(kmid,kmid,e) ,1:2); 
      if ( midpoint(1)<xmin | midpoint(1)>xmax | midpoint(2)<ymin | midpoint(2)>ymax ), continue; end
      indx(ip+1:ip+blocksize,:) = indx0 + (e-1)*NGLL^2;
      ip = ip + blocksize;
    end
    if ip<size(indx,1), indx = indx(1:ip,4); end

  else
    for e=1:NEL,
      ip = (e-1)*blocksize;
      indx(ip+1:ip+blocksize,:) = indx0 + (e-1)*NGLL^2;
    end
  end

end
