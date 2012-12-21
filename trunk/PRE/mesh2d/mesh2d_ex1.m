% MESH2D_EX1	Mesh for a shallow layer over half-space with dipping fault
%
% Shallow low velocity layer, thickness = 5 km
% The fault intersects the free surface at (0,0)
% It has two segments:
%   Shallow segment :  dip 70 degrees   length 8 km
%   Deep segment    :  dip 50 degrees   length 30 km
%

%-- you can modify the following parameters:

LAYER_DEPTH = 5e3; 		% thickness of shallow layer
BOX = [-30e3 50e3 -50e3 0]; 	% domain limits [xmin xmax zmin zmax]
h = 500;			% element size along fault (approximate)
%h = 1000;			% element size along fault (approximate)
%h = 2000;			% element size along fault (approximate)

% fault geometry parameters
dip1 = 70*pi/180;
dip3 = 50*pi/180;
L1 = 8e3;
L2 = 30e3;
DYNFLT = 1;	% 1 = dynamic fault (split nodes), 0 = fault is welded (no split nodes)

%-- you dont't need to modify anything below this line

% points on the fault
p1 = [0 0];
p3 = p1 + [cos(dip1) -sin(dip1)]*L1;
p2 = [p3(1)+(p1(1)-p3(1))*(-LAYER_DEPTH-p3(2))/(p1(2)-p3(2)), -LAYER_DEPTH]; % fault/layer intersection
p4 = p3 + [cos(dip3) -sin(dip3)]*L2;
p5 = [p4(1) BOX(3)];
fault = [p5; p4; p3; p2; p1];

% points on the left boundary
left = [[BOX(1) BOX(3)]; ...
       [BOX(1) p4(2)]; ...
       [BOX(1) p3(2)]; ...
       [BOX(1) -LAYER_DEPTH]; ...
       [BOX(1) BOX(4)]];

% points on the right boundary
right = [[BOX(2) BOX(3)]; ...
       [BOX(2) p4(2)]; ...
       [BOX(2) p3(2)]; ...
       [BOX(2) -LAYER_DEPTH]; ...
       [BOX(2) BOX(4)]];

figure(1)
clf
plot( fault(:,1),fault(:,2),'o', left(:,1),left(:,2),'o', right(:,1),right(:,2),'o' )
axis equal
axis(BOX)

% number of elements along x :
% such that horizontal element size at mid fault depth = h
NELX(1) = ceil(abs((p1(1)+p4(1))/2-BOX(1))/h);	% in the left side
NELX(2) = ceil(abs(BOX(2)-(p1(1)+p4(1))/2)/h);	% in the right side

% number of elements along z :
NELZ(4) = ceil(LAYER_DEPTH/sin(dip1)/h);	% in the shallow layer
NELZ(3) = ceil((L1-LAYER_DEPTH/sin(dip1))/h);	% from bottom of shallow layer to fault kink
NELZ(2) = ceil(L2/h);				% from fault kink to fault end
NELZ(1) = ceil((p4(2)-p5(2))/h);		% from fault end to bottom of domain


%-- mesh the block on the left side of the fault

% Define the 4 boundaries:
curves{1} = sample_segments([[BOX(1) BOX(3)]; p5],NELX(1)); 	% bottom 
curves{2} = sample_segments(fault,NELZ);		 	% right 
curves{3} = sample_segments([[BOX(1) BOX(4)]; p1],NELX(1)); 	% top 
curves{4} = sample_segments(left,NELZ); 			% left 

% Generate the mesh
domain(1) = mesh2d_tfi(curves,4);

% tag the shallow layer
e = [sum(NELZ(1:3))*NELX(1)+1 : sum(NELZ)*NELX(1)];
domain(1).etag(e) = 2;

% tag the elements close to the fault
if DYNFLT
  e = [ NELZ(1) : sum(NELZ(1:3)) ] *NELX(1);
  domain(1).etag(e) = 3;
  e = [ sum(NELZ(1:3))+1 : sum(NELZ) ] *NELX(1);
  domain(1).etag(e) = 4;
end

%-- mesh the block on the right side of the fault

% Define the 4 boundaries:
curves{4} = curves{2}; 						% left 
curves{1} = sample_segments([p5; [BOX(2) BOX(3)]],NELX(2)); 	% bottom 
curves{2} = sample_segments(right,NELZ); 			% right 
curves{3} = sample_segments([p1; [BOX(2) BOX(4)]],NELX(2)); 	% top 

% Generate the mesh
domain(2) = mesh2d_tfi(curves,4);

% tag the shallow layer
e = [sum(NELZ(1:3))*NELX(2)+1 : sum(NELZ)*NELX(2)];
domain(2).etag(e) = 2;

% tag the elements close to the fault
if DYNFLT
  e = 1+ [ NELZ(1)-1 : sum(NELZ(1:3))-1 ] *NELX(2);
  domain(2).etag(e) = 3;
  e = 1+ [ sum(NELZ(1:3)) : sum(NELZ)-1 ] *NELX(2);
  domain(2).etag(e) = 4;
end

%-- plot the two meshes
mesh2d_plot(domain,1)

%-- merge the two domains into a single mesh
if DYNFLT
  mtab = [ [-1; -5; -3; -4], [-1; -2; -3; -6] ];
else
  mtab = [ [-1; 2; -3; -4], [-1; -2; -3; 1] ];
end
mesh = mesh2d_merge(domain,mtab)
figure(2)
clf
mesh2d_plot(mesh,2)

%-- export the mesh file
mesh2d_write(mesh,'ex1');
