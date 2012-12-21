% SEM2D_READ_SPECGRID reads a spectral element grid
%
% SYNTAX	grid = sem2d_read_specgrid(datadir)
%
% INPUT		datadir	data directory name, default is current directory
%
% OUTPUT	grid 	spectral element grid structure containing:
%			nelem	number of elements
%			ngnod	number of control nodes per element (4 or 9)
%			npgeo	total number of control nodes
%			coorg	coordinates of control nodes
%			knods	connectivity table of the finite element mesh, 
%				size=[nelem,ngnod]
%				knods(k,e) is the index in the list of control nodes
%				of the k-th control node of the e-th element
%			npoin	number of global GLL nodes
%			ngll	number of GLL nodes per element side
%			coord	coordinates of global GLL nodes, size (npoin,2)
%			ibool	local to global numbering, size=[ngll,ngll,nelem]
% 				needed to lookup node coordinates for fields stored 
%				in element-wise format.
% 				ibool(i,j,e) is the global index of the (i,j)-th 
%				local GLL node of the e-th element. 
%			x	coordinates of GLL points in the reference 
%				segment [-1:1], size=ngll
%               	w    	GLL quadrature weights, size=ngll
%               	h(:,:)  derivatives of Lagrange polynomials at the 
%				GLL nodes: h(i,j) = h'_i( x(j) )
%
% NOTE	The spectral element grid is built from an initial finite element mesh
% 	by subdividing each element into a subgrid of ngll*ngll internal nodes,
% 	the Gauss-Lobatto-Legendre (GLL) nodes.
% 	The local subgrids are structured, possibly deformed,
% 	but irregularly spaced (GLL nodes cluster at the element edges). 
%
function g = sem2d_read_specgrid(datadir)

if exist('datadir','var') 
  if datadir(end) ~= '/', datadir = [datadir '/']; end
else
  datadir='';
end

% read mesh information
[g.nelem,g.npgeo,g.ngnod,g.npoin,g.ngll] = ...
  textread([datadir 'grid_sem2d.hdr'],'%u%u%u%u%u','headerlines',1);

% read finite element grid
g.coorg = load('-ascii',[datadir 'MeshNodesCoord_sem2d.tab']);
g.knods = load('-ascii',[datadir 'ElmtNodes_sem2d.tab']);

% read spectral element grid
fid=fopen([datadir 'coord_sem2d.dat']);
g.coord = fread(fid,[2,g.npoin],'single')' ; 
fclose(fid);

fid=fopen([datadir 'ibool_sem2d.dat']); 
g.ibool = fread(fid,inf,'int'); 
fclose(fid);
g.ibool = reshape(g.ibool,[g.ngll,g.ngll,g.nelem]);

% read GLL information
fid=fopen([datadir 'gll_sem2d.tab']);
data=fscanf(fid,'%f',[g.ngll,g.ngll+2]);
fclose(fid);
g.x=data(:,1);		% GLL nodes in reference segment [-1:1]
g.w=data(:,2);		% GLL quadrature weights
g.h=data(:,3:end)';	% Derivatives of Lagrange polynomials at the GLL nodes
			%  h(i,j) = h'_i( x(j) )
