% MESH2D_ROTATE rotates the node coordinates of a 2D mesh
%
% SYNTAX	rmesh = mesh2d_rotate(mesh,theta)
%
% INPUTS	mesh	a mesh2d structure (see MESH2D_TFI)
%		theta 	rotation angle, clockwise, in radians
%
% OUTPUTS	rmesh	new mesh2d structure with rotated nodes
%
function mesh2 = mesh2d_rotate(mesh1,theta)
mesh2 = mesh1;
c = cos(theta);
s = sin(theta);
mesh2.coor(1,:) = c*mesh1.coor(1,:) - s*mesh1.coor(2,:);
mesh2.coor(2,:) = s*mesh1.coor(1,:) + c*mesh1.coor(2,:);
