% SEM2D_READ_SPECGRID reads a spectral element grid
%
% SYNTAX	grid = sem2d_read_specgrid(datadir)
%
% INPUT		datadir	data directory name, default is current directory
%
% OUTPUT	grid 	spectral element grid structure containing:
%			nelem	number of elements
%			npoin	number of global GLL nodes
%			ngll	number of GLL nodes per element side
%			coord	coordinates of global GLL nodes, size (npoin,2)
%			ibool	local to global numbering, size=[ngll,ngll,nelem]
% 				needed to lookup node coordinates for fields stored 
%				in element-wise format.
% 				ibool(i,j,e) is the global index of the (i,j)-th 
%				local GLL node of the e-th element. 
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

[g.nelem,g.npoin,g.ngll] = textread([datadir 'grid_sem2d.hdr'],'%u%*u%*u%u%u','headerlines',1);

fid=fopen([datadir 'coord_sem2d.dat']);
g.coord = fread(fid,[2,g.npoin],'single')' ; 
fclose(fid);

fid=fopen([datadir 'ibool_sem2d.dat']); 
g.ibool = fread(fid,inf,'int'); 
fclose(fid);
g.ibool = reshape(g.ibool,[g.ngll,g.ngll,g.nelem]);

