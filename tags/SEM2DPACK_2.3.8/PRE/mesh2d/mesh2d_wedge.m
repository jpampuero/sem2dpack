% MESH2D_WEDGE  Generates a mesh for a triangular wedge domain
%		by merging three structured meshes
%
% SYNTAX	mesh = mesh2d_wedge(X,Y,NEL,Q)
%               mesh = mesh2d_wedge(curves,Q)
%
% INPUTS	X,Y	position of the vertices of the triangle, length=3
%		NEL	half the number of elements per triangle edge
% 		curves 	cell array of length 3 containing the nodes along the
%			3 boundaries of the domain: curves{k}(i,1:2) is the
%			location (x,y) of the i-th node along the k-th boundary.
%			Boundaries and nodes are numbered counter-clockwise.
%			All boundaries must have the same odd number of nodes.
%			The last node of a boundary must coincide with the first
%			node of the next boundary.
%		Q	number of nodes per quad element: 4 or 9 (default=4)
%
% OUTPUTS	mesh	a mesh2d structure (see MESH2D_TFI)
%
function mesh = mesh2d_wedge(varargin)

if iscell(varargin{1})
  mesh = mesh2d_wedge_tfi(varargin{:});
else
  mesh = mesh2d_wedge_quad(varargin{:});
end


%-----------
%function mesh = mesh2d_wedge_quad(X,Y,NEL,Q)
function mesh = mesh2d_wedge_quad(varargin)

X = varargin{1};
Y = varargin{2};
NEL = varargin{3};
if nargin==4, Q=varargin{4}; end

if length(X)~=3 | length(Y)~=3, error('X and Y must have length = 3'); end
if ~exist('Q','var'), Q=4; end

X(4) = (X(1)+X(2))/2;
X(5) = (X(2)+X(3))/2;
X(6) = (X(3)+X(1))/2;
X(7) = (X(1)+X(2)+X(3))/3;

Y(4) = (Y(1)+Y(2))/2;
Y(5) = (Y(2)+Y(3))/2;
Y(6) = (Y(3)+Y(1))/2;
Y(7) = (Y(1)+Y(2)+Y(3))/3;

pts = [1 4 7 6];
domain(1) = mesh2d_quad(X(pts),Y(pts),NEL,NEL,Q);

pts = [4 2 5 7];
domain(2) = mesh2d_quad(X(pts),Y(pts),NEL,NEL,Q);

pts = [6 7 5 3];
domain(3) = mesh2d_quad(X(pts),Y(pts),NEL,NEL,Q);

mtab = [ [-1; 2; 3; -3], [-1; -2; 3; 1], [1; 2; -2; -3] ];
mesh = mesh2d_merge(domain,mtab);

%--------
% function mesh = mesh2d_wedge_tfi(curves,Q)
function mesh = mesh2d_wedge_tfi(varargin)

% assumes curves are numbered counterclockwise
% and nodes in each curve are ordered counterclockwise

curves = varargin{1};
if nargin>1, Q=varargin{2}; end
if nargin>2, nel=varargin{3}; end
if nargin>3, P7=varargin{4}; end

if ~exist('Q','var'), Q=4; end

N1 = size(curves{1},1)-1;
N2 = size(curves{2},1)-1;
N3 = size(curves{3},1)-1;

if ~exist('nel','var') 
  nel = floor((N1-1)/2);
end

nels(1) = nel;
nels(2) = N1-nels(1);
nels(3) = N2-nels(1);

if N3~=nels(2)+nels(3)
  N3_is = N3
  N3_should_be = nels(2)+nels(3)
  error('Incompatible boundaries'); 
end

P1 = curves{1}(1,:);
P2 = curves{2}(1,:);
P3 = curves{3}(1,:);
P4 = curves{1}(nels(1)+1,:);
P5 = curves{2}(nels(3)+1,:);
P6 = curves{3}(nels(2)+1,:);
if ~exist('P7','var'), P7 = (P4+P5+P6)/3; end

segm1 = sample_segments([P4;P7],nels(3));
segm2 = sample_segments([P7;P5],nels(2));
segm3 = sample_segments([P6;P7],nels(1));

cur{1} = curves{1}(1:nels(1)+1,:);
cur{2} = segm1;
cur{3} = segm3;
cur{4} = curves{3}(end:-1:nels(2)+1,:);
domain(1) = mesh2d_tfi(cur,Q);

cur{1} = curves{1}(nels(1)+1:end,:);
cur{2} = curves{2}(1:nels(3)+1,:);
cur{3} = segm2;
cur{4} = segm1;
domain(2) = mesh2d_tfi(cur,Q);

cur{1} = segm3;
cur{2} = segm2;
cur{3} = curves{2}(end:-1:nels(3)+1,:);
cur{4} = curves{3}(nels(2)+1:-1:1,:);
domain(3) = mesh2d_tfi(cur,Q);

mtab = [ [-1; 2; 3; -3], [-1; -2; 3; 1], [1; 2; -2; -3] ];
mesh = mesh2d_merge(domain,mtab);

