% MESH2D_QUAD	generates a structured mesh for a quadrilateral domain
%	
% SYNTAX	To mesh a general quadrilateral domain:
%		  mesh = mesh2d_quad(X,Y,NELX,NELY,Q)
%		To mesh a rectangular domain:
% 		  mesh = mesh2d_quad(LX,LY,NELX,NELY,Q)
%
% INPUTS	X(1:4)	x-coordinates of domain vertices
%		Y(1:4) 	y-coordinates of domain vertices
%			The numbering convention for vertices is:
%			
%            		  4---------3
%               	  |         |
%               	  |         |
%            		  1---------2
%
% 		LX	size of rectangular domain in the x direction
% 		LY	size of rectangular domain in the y direction
%			The rectangular domain is
%
%            		(0,LY)-----(LX,LY)
%               	  |           |
%               	  |           |
%            		(0,0)------(LX,0)
%			
%		NELX	number of elements in the direction of boundaries #1 and #3
%		NELY	number of elements in the direction of boundaries #2 and #4 
%			The numbering convention for boundaries is
%			
%            		  .----3----.
%               	  |         |
%             		  4         2
%               	  |         |
%            		  '----1----'
%
%		Q	number of nodes per element: 4 or 9
% 		
% OUTPUTS	mesh	a mesh2d structure (see MESH2D_TFI)
%

function mesh = mesh2d_quad(X,Y,NELX,NELY,Q)

if ~((length(X)==1 & length(Y)==1) | (length(X)==4 & length(Y)==4))
  error('The lengths of X and Y should be both 1 or 4');
elseif length(X)==1
  X = [0 X X 0];
  Y = [0 0 Y Y];
end

if ~exist('Q','var'), Q=4; end
if Q==4
  nex= NELX;
  ney= NELY;
elseif Q==9
  nex= 2*NELX;
  ney= 2*NELY;
else
  error('Only Q4 and Q9 elements are implemented')
end

curves{1} = sample_segments([ [X(1) Y(1)]; [X(2) Y(2)] ],nex);
curves{2} = sample_segments([ [X(2) Y(2)]; [X(3) Y(3)] ],ney);
curves{3} = sample_segments([ [X(4) Y(4)]; [X(3) Y(3)] ],nex);
curves{4} = sample_segments([ [X(1) Y(1)]; [X(4) Y(4)] ],ney);

mesh = mesh2d_tfi(curves,Q);
