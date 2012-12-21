% MESH2D_TFI	generates a structured mesh by transfinite interpolation
% 		for a deformed quadrilateral domain
%
% SYNTAX	mesh = mesh2d_tfi(curves,Q)
% 		mesh = mesh2d_tfi('demo')
%
% INPUT		curves 	cell array of length 4 containing the nodes along the
%			4 boundaries of the domain: curves{k}(i,1:2) is the
%			location (x,y) of the i-th node along the k-th boundary.
%			Boundaries must be numbered counter-clockwise.
%			Facing boundaries must have the same number of nodes.
%			Nodes must be ordered counter-clockwise in boundaries 1 and 2,
%			clockwise in boundaries 3 and 4.
%			The last node of a boundary must coincide with the first
%			node of the next boundary.
%		Q	quad element type, by number of nodes: 4 or 9
%			If Q=9, the number of nodes along each boundary must be odd.
%
% OUTPUT	mesh	mesh structure, contains:
%			coor	size = [2, number of mesh nodes]
%				coordinates (x,y) of the mesh nodes
%			enod	size = [Q, number of elements]
%				enod(i,e) is the global index of the i-th node
%				of the e-th element.
%			etag	length = number of elements
%				etag(e) is the material tag of the e-th element
%			bnds	cell array of length = number of boundaries
%				The eb-th boundary element of the k-th boundary
%				is the bnds{k}(2,eb)-th edge of the bnds{k}(1,eb)-th 
%				element of mesh.
%
%				Inside each element the nodes and edges (numbers in [])
%				are locally numbered as
%			
%                                   [3]
%		            	4----7----3
%               		|         |
%            		  [4]   8    9    6   [2]
%           		    	|         |
%           		 	1----5----2
%                                   [1]
%
% EXAMPLE	m=mesh2d_tfi('demo');
%
function mesh = mesh2d_tfi(curves,Q)

% check inputs
if isstr(curves) && strcmp('demo',curves)
  [curves,Q] = demo_curves();
  DEMO=1;
else
  DEMO=0;
end
if ~iscell(curves), error('Input curves must a cell array of length 4'), end
if length(curves)~=4, error('Input curves must a cell array of length 4'), end
if size(curves{1},1)~=size(curves{3},1), error('curves{1} and curves{3} must have same size'), end
if size(curves{2},1)~=size(curves{4},1), error('curves{2} and curves{4} must have same size'), end
if ~exist('Q','var'), Q=4; end
if Q~=4 & Q~=9,  error('Only Q4 and Q9 elements are implemented'), end

npx = size(curves{1},1);
npy = size(curves{2},1);
if Q==4
  NELX = npx-1;
  NELY = npy-1;
else
  NELX = (npx-1)/2;
  NELY = (npy-1)/2;
  if NELX~=floor(NELX) | NELY~=floor(NELY) , error('Number of nodes on curves must be odd if Q=9'), end
end

mesh.coor = tfi_COOR(curves);
mesh.enod = structured_ENOD(NELX,NELY,Q); 	% compute element table
mesh.etag = ones(1,NELX*NELY);
mesh.bnds = structured_BNDS(NELX,NELY); 	% compute boundary table

if DEMO, mesh2d_plot(mesh); end

%--------------------------------------------
function [curves,Q] = demo_curves()

npx = 13;
npy = 7;
xi = linspace(0,1,npx)';
eta = linspace(0,1,npy)';
curves{1} = [xi, -0.1*sin(pi*xi).^2];
curves{2} = [1+0.1*sin(pi*eta).^2, eta];
curves{3} = [xi, 1+0.2*sin(pi*xi).^2];
curves{4} = [-0.2*sin(pi*eta).^2, eta];
Q=9;

%--------------------------------------------
% compute coordinates by transfinite interpolation
function coor = tfi_COOR(curves)

npx = size(curves{1},1);
npy = size(curves{2},1);

xi = linspace(0,1,npx)';
eta= linspace(0,1,npy);

shape1 = [1-xi , xi ];
shape2 = [1-eta; eta];

for k=1:2,
  xy = shape1 * [curves{4}(:,k)'; curves{2}(:,k)'] ...
     + [curves{1}(:,k), curves{3}(:,k)] * shape2 ...
     - shape1 * [curves{1}(1,k)   curves{3}(1,k); ...
                 curves{1}(end,k) curves{3}(end,k)] * shape2;
  coor(k,:) = reshape(xy,1,npx*npy);
end


%--------------------------------------------
% compute element table
% enod 	size=[Q,NELX*NELY]
%
function enod = structured_ENOD(NELX,NELY,Q)

if ~exist('Q','var'), Q=4; end
if Q==4
  npx= NELX+1;
  npy= NELY+1;
 % node coordinates in the reference element
  irel = [ 0 1 1 0 ]';
  jrel = [ 0 0 1 1 ]';
elseif Q==9
  npx= 2*NELX+1;
  npy= 2*NELY+1;
  irel = [ -1  1 1 -1  0 1 0 -1 0 ]';
  jrel = [ -1 -1 1  1 -1 0 1  0 0 ]';
else
  error('Only Q4 and Q9 elements are implemented')
end

% elements numbered from left to right, from bottom to top
ie = repmat([1:NELX],1,NELY);
je = repmat([1:NELY],NELX,1);
je=je(:)';
if Q==9
  ie = 2*ie;
  je = 2*je;
end
ip = repmat(ie,Q,1) + repmat(irel,1,NELX*NELY);
jp = repmat(je,Q,1) + repmat(jrel,1,NELX*NELY);
enod = sub2ind([npx,npy],ip,jp); 	% size [Q,NELX*NELY]

%--------------------------------------------
% compute boundary tables
% bnds	cell array, length=4 (number of boundaries)
%	ordered as: 1=bottom, 2=right, 3=top, 4=left
%       size of bnds{k} = [2, number of elements in k-th boundary]
%	bnds{k}(1,:) = element id
%	bnds{k}(2,:) = face id, in local frame: 1=bottom, 2=right, 3=top, 4=left
%
function bnds = structured_BNDS(NELX,NELY)

bnds{1} = [ (1:NELX) ; ones(1,NELX) ];
bnds{2} = [ (1:NELY)*NELX ; repmat(2,[1,NELY]) ];
bnds{3} = [ (NELX*(NELY-1)+1:NELY*NELX); repmat(3,[1,NELX]) ];
bnds{4} = [ (0:NELY-1)*NELX+1 ; repmat(4,[1,NELY]) ];

