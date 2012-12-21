% MESH2D_TRANSLATE translates the node coordinates of a 2D mesh
%
% SYNTAX	tmesh = mesh2d_rotate(mesh,t)
%
% INPUTS	mesh	a mesh2d structure (see MESH2D_TFI)
%		t(1:2)  translation vector
%
% OUTPUTS	tmesh	new mesh2d structure with translated nodes
%

function mesh2 = mesh2d_translate(mesh1,t)
mesh2 = mesh1;
mesh2.coor(1,:) = mesh1.coor(1,:)+t(1);
mesh2.coor(2,:) = mesh1.coor(2,:)+t(2);
