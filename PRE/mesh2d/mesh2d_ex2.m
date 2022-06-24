% MESH2D_EX2	Two horizontal faults cutting through a rectangular domain
%
% Boundary tags:
% Bottom fault: 5 = bottom side, 6 = top side
% Top fault:    7 = bottom side, 8 = top side
%
% Element tags:
% Elements touching the faults: tag = 2
% Other elements: tag = 1

%-- you can modify the following parameters:
BOX = [-1 1 -1 1]*2000; % domain limits [xmin xmax zmin zmax]
Z_F1 = -40; 	% vertical position of bottom fault
Z_F2 = -Z_F1; 	% vertical position of top fault
h = 10;			% element size

%-- no need to modify anything below this line

% make sure fault 1 is the bottom one
% swap them if needed
if (Z_F2 < Z_F1)
  tmp = Z_F1;
  Z_F1 = Z_F2;
  Z_F2 = tmp;
end

% box vertices
p1 = [BOX(1) BOX(3)]; % bottom left
p2 = [BOX(2) BOX(3)]; % bottom right
p3 = [BOX(2) BOX(4)]; % top right
p4 = [BOX(1) BOX(4)]; % top left

% end points of fault 1
pF1_1 = [BOX(1) Z_F1]; % F1 left
pF1_2 = [BOX(2) Z_F1]; % F1 right

% end points of fault 2
pF2_1 = [BOX(1) Z_F2]; % F2 left
pF2_2 = [BOX(2) Z_F2]; % F2 right

% three sub-domains separated by the two faults:

% 4 -------------------------- 3
% |                            |
% F2_1 -------------------- F2_2
% |                            |
% F1_1 -------------------- F1_2
% |                            |
% 1 -------------------------- 2

NELX = ceil((BOX(2)-BOX(1))/h);

% bottom sub-domain
ps = [p1; p2; pF1_2; pF1_1];
NELZ_1 = ceil((Z_F1-BOX(3))/h);
domain(1) = mesh2d_quad(ps(:,1),ps(:,2),NELX,NELZ_1);

% middle sub-domain
ps = [pF1_1; pF1_2; pF2_2; pF2_1];
NELZ_2 = ceil((Z_F2-Z_F1)/h);
domain(2) = mesh2d_quad(ps(:,1),ps(:,2),NELX,NELZ_2);

% top sub-domain
ps = [pF2_1; pF2_2; p3; p4];
NELZ_3 = ceil((BOX(4)-Z_F2)/h);
domain(3) = mesh2d_quad(ps(:,1),ps(:,2),NELX,NELZ_3);

% tag the elements touching the faults
% to create Kelvin-Voigt viscous layers
e = sub2ind([NELX NELZ_1], [1:NELX]',NELZ_1*ones(NELX,1));
domain(1).etag(e) = 2;
e = sub2ind([NELX NELZ_2], [1:NELX]',ones(NELX,1));
domain(2).etag(e) = 2;
e = sub2ind([NELX NELZ_2], [1:NELX]',NELZ_2*ones(NELX,1));
domain(2).etag(e) = 2;
e = sub2ind([NELX NELZ_3], [1:NELX]',ones(NELX,1));
domain(3).etag(e) = 2;

% merge domains into a single mesh
mtab = [ [ -1 -6 -8 ] ; ...
         [ -2 -2 -2 ] ; ...
         [ -5 -7 -3 ] ; ...
         [ -4 -4 -4 ] ];
mesh = mesh2d_merge(domain,mtab)

figure(2)
clf
mesh2d_plot(mesh,2)
axis equal

%-- export the mesh file
mesh2d_write(mesh,'ex2');
