% MESH2D_EX0	2D mesh for a vertical cross-section across a vertical strike-slip fault

%-- you can modify the following parameters:
FAULT_BOTTOM = 15e3; 	% depth of the bottom tip of the fault
FAULT_X = 0;
BOX = [-15e3 15e3 -30e3 0]; 	% domain limits [xmin xmax zmin zmax]
h = 500;			% element size

%-- you don't need to modify anything below this line

NELX_L = ceil((FAULT_X-BOX(1))/h);
NELX_R = ceil((BOX(2)-FAULT_X)/h);
NELZ_T = ceil((BOX(4)+FAULT_BOTTOM)/h);
NELZ_B = ceil((-FAULT_BOTTOM-BOX(3))/h);

% box vertices
p1 = [BOX(1) BOX(3)];
p2 = [BOX(2) BOX(3)];
p3 = [BOX(2) BOX(4)];
p4 = [BOX(1) BOX(4)];

% points on the fault
p7 = [0 0];
p9 = [0 -FAULT_BOTTOM];

% points on box edges, projecting fault bottom tip
p5 = [0 BOX(3)];
p6 = [BOX(2) -FAULT_BOTTOM];
p8 = [BOX(1) -FAULT_BOTTOM];

% bottom left domain
ps = [p1; p5; p9; p8];
domain(1) = mesh2d_quad(ps(:,1),ps(:,2),NELX_L,NELZ_B);

% bottom right domain
ps = [p5; p2; p6; p9];
domain(2) = mesh2d_quad(ps(:,1),ps(:,2),NELX_R,NELZ_B);

% top right domain
ps = [p9; p6; p3; p7];
domain(3) = mesh2d_quad(ps(:,1),ps(:,2),NELX_R,NELZ_T);

% top left domain
ps = [p8; p9; p7; p4];
domain(4) = mesh2d_quad(ps(:,1),ps(:,2),NELX_L,NELZ_T);

% tag the elements close to the fault
% to create a Kelvin-Voigt viscous layer
e = sub2ind([NELX_R NELZ_B], NELX_L,NELZ_B);
domain(1).etag(e) = 2;
e = sub2ind([NELX_R NELZ_B], 1,NELZ_B);
domain(2).etag(e) = 2;
e = sub2ind([NELX_R NELZ_T], ones(NELZ_T,1),[1:NELZ_T]');
domain(3).etag(e) = 2;
e = sub2ind([NELX_L NELZ_T], NELX_L*ones(NELZ_T,1),[1:NELZ_T]');
domain(4).etag(e) = 2;

% merge domains into a single mesh
mtab = [ [ -1 -1  2  1] ; ...
         [  2 -2 -2 -5] ; ...
         [  4  3 -3 -3] ; ...
         [ -4  1 -6 -4] ];
mesh = mesh2d_merge(domain,mtab)
figure(2)
clf
mesh2d_plot(mesh,2)
axis equal

%-- export the mesh file
mesh2d_write(mesh,'ex0');
