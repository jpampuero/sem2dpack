% MESH2D_CIRC_HOLE generates a mesh for a square domain with a central circular hole 
%
% SYNTAX	mesh = mesh2d_circ_hole(L,R,NPL,NPR)
% 		mesh = mesh2d_circ_hole('demo')
%
% INPUTS	L	side length of the square domain
%		R	radius of the inner hole
%		NPL	number of elements per side
%		NPR	number of elements from side to hole
%
% OUTPUTS	mesh	a mesh structure
%
% EXAMPLE	m=mesh2d_circ_hole('demo')
%
function mesh = mesh2d_circ_hole(L,R,NPL,NPR)

if isstr(L) && strcmp('demo',L)
  DEMO = 1;
  L = 10;
  R = 3;
  NPL = 6;
  NPR = 3;
else 
  DEMO =0;
end

%-- 1. Generate four structured domains, each one 
% connects a 1/4-arc of the central hole to one exterior face.
% Define each domain separately: 
%        coor = coordinates of each node
%        enod = index of nodes of each element
%        bnds = bulk element and face for each boundary element

% define domain#1
npl=2*NPL+1; % number of nodes, Q9 elements
curve{1} = sample_segments([[-L/2 -L/2];[L/2 -L/2]],npl-1);
theta = linspace(-3*pi/4,-pi/4,npl)';
curve{3} = [ R*cos(theta) R*sin(theta) ];
curve{2} = sample_segments([curve{1}(npl,:); curve{3}(npl,:)],2*NPR);
curve{4} = sample_segments([curve{1}(1,:); curve{3}(1,:)],2*NPR);
domain(1) = mesh2d_tfi(curve,9);

% other domains are obtained by rotation
% must preserve orientation of the elements
domain(2)= mesh2d_rotate(domain(1),pi/2);
domain(3)= mesh2d_rotate(domain(1),pi);
domain(4)= mesh2d_rotate(domain(1),3*pi/2);

if DEMO
  mesh2d_plot(domain); 
end

%-- 2. Merge
MergeTo(:,1) = [ -1; 2;-5; 4 ]; 
MergeTo(:,2) = [ -2; 3;-5; 1 ]; 
MergeTo(:,3) = [ -3; 4;-5; 2 ]; 
MergeTo(:,4) = [ -4; 1;-5; 3 ]; 
mesh = mesh2d_merge(domain,MergeTo);

