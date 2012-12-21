% MESH2D_WEDGE  Generates a mesh for a triangular wedge domain
%		by merging three structured meshes
%
% SYNTAX	mesh = mesh2d_wedge(X,Y,NEL,Q)
%
% INPUTS	X,Y	position of the vertices of the triangle, length=3
%		NEL	half the number of elements per triangle edge
%		Q	number of nodes per quad element: 4 or 9 (default=4)
%
% OUTPUTS	mesh	a mesh2d structure (see MESH2D_TFI)
%
function mesh = mesh2d_wedge(X,Y,NEL,Q)

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
mesh = mesh2d_merge(domain,mtab)
