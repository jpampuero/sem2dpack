% This matlab script plots snapshots from SEM2D output.
% You must set FNAME and ISNAP below.

% select a field : 
%			PSV			SH
% 	displacement	dx,dz			dy
%	velocity	vx,vz			vy
%	acceleration	ax,az			ay
%	strain		e11,e22,e12 		e13,e23
%	stress		s11,s22,s33,s12		s13,s23
FNAME = 'vy';
% select snapshot indices:
ISNAP = [1:20];

%------------------------------------------------------

% The spectral element grid is built from an initial "macro-mesh"
% (essentially a finite element mesh, generated for instance by EMC2)
% by subdividing each element into a subgrid of ngll*ngll internal nodes,
% the Gauss-Lobatto-Legendre (GLL) nodes.
% These local subgrids are structured (deformed cartesian grids),
% but not regularly spaced (GLL nodes cluster at the element edges). 

%-- Prepare grid data :

% Grid parameters :
%	nelem	number of elements
%	npgeo	total number of points of the macro-mesh	
%	ngnod	number of nodes per element in the macro-mesh
% 	npoin	total number of GLL nodes
% 	ngll	number of GLL nodes per element edge
[nelem,npgeo,ngnod,npoin,ngll] = ...
  textread('grid_sem2d.hdr','%u%u%u%u%u','headerlines',1);

% 'd', 'v' and 'a' are stored in "global" format, node-by-node.
% 'e' and 's' are stored in "local" format, element-by-element.
% For instance, vx(:) is a vector of length npoin
% and e11(:,:,:) is a 3d order matrix of dimensions ngll*ngll*nelem.
IsNodal = isempty(findstr(FNAME(1),'es'));

% coord(npoin,2) = coordinates of the GLL nodes
fid=fopen('coord_sem2d.dat'); coord = fread(fid,[2,npoin],'single')' ; fclose(fid);

% ibool(ngll,ngll,nelem) = local to global numbering
% This table is needed to lookup node coordinates for fields stored in "local" format.
% ibool(i,j,e) is the global index of the (i,j) GLL node local to element e. 
fid=fopen('ibool_sem2d.dat'); ibool = fread(fid,inf,'int'); fclose(fid);
ibool = reshape(ibool,[ngll,ngll,nelem]);


%-- Initialize data for snapshots:
%   build the cell data for Matlab's PATCH function
if IsNodal
  indx = Init2dSnapshot(ibool);
else
  [indx,coord] = Init2dSnapshot(ibool,coord);
end

%-- Plot 
for is=ISNAP,
  fname = sprintf('%s_%3.3u_sem2d.dat',FNAME,is);
  if ~exist(fname,'file'), break, end
  fid = fopen(fname); field = fread(fid,inf,'single'); fclose(fid);
  if ~IsNodal, field=reshape(field,[ngll,ngll,nelem]); end
  clf
  Plot2dSnapshot(field,coord,indx);
  pause
end
