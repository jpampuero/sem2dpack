% build mesh for rate state cycle simulation
%
% Two-sided faults
%
% a finite single fault in the middle of the domain
%
%   | --------------------------------------------------|
%   |                                                               |
%   |                                                               |
%   |                     ----------------                      |
%   |                                                               |
%   |                                                               |
%   | --------------------------------------------------|
%

addpath(genpath('../../PRE/'));

close all;

%
Lx    =   91;
Lz    =   120;

%-- you can modify the following parameters:
FAULT_L    = 45; % fault length 	
FAULT_Z    = 0;   % fault z location 
FZW           = 0;   % fault weak zone thickness each side
KVW           = 0;  % width of KV material on each side
tag_fz         = 2;   % tag for fault zone
tag_kv        = 3;   % tag for KV damper, tag_kv=0, no KV damper
h                 = 1;    % element size

% ----------------------------- Do not modify --------------------------
% create vertices
p1   = [-Lx/2,             -Lz/2];
p2   = [-FAULT_L/2,  -Lz/2];
p3   = [FAULT_L/2,   -Lz/2];
p4   = [Lx/2,             -Lz/2];

p5   = [-Lx/2,              0];
p6   = [-FAULT_L/2,   0];
p7   = [FAULT_L/2,    0];
p8   = [ Lx/2,               0];

p9   = [-Lx/2,              Lz/2];
p10 = [-FAULT_L/2,   Lz/2];
p11 = [FAULT_L/2,    Lz/2];
p12 = [Lx/2,               Lz/2];

% number of elements for left, middle, right columns
NELX_L  =  ceil((Lx - FAULT_L)/2/h);
NELX_M = ceil(FAULT_L/h);
NELX_R =  NELX_L;

NELZ_T = ceil(Lz/2/h);
NELZ_B = NELZ_T;

% create 6 domains

% LEFT-BOTTOM
ps = [p1; p2; p6; p5];
domain(1) = mesh2d_quad(ps(:,1), ps(:,2), NELX_L, NELZ_B);

% MID-BOTTOM
ps = [p2; p3; p7; p6];
domain(2) = mesh2d_quad(ps(:,1), ps(:,2), NELX_M, NELZ_B);

% RIGHT-BOTTOM
ps = [p3; p4; p8; p7];
domain(3) = mesh2d_quad(ps(:,1), ps(:,2), NELX_R, NELZ_B);

% RIGHT-TOP
ps = [p7; p8; p12; p11];
domain(4) = mesh2d_quad(ps(:,1), ps(:,2), NELX_R, NELZ_T);

% MID-TOP
ps = [p6; p7; p11; p10];
domain(5) = mesh2d_quad(ps(:,1), ps(:,2), NELX_M, NELZ_T);

% LEFT-TOP
ps = [p5; p6; p10; p9];
domain(6) = mesh2d_quad(ps(:,1), ps(:,2), NELX_L, NELZ_T);

% tag the fault zone
if FZW >= h
    iFZW = ceil(FZW/h);
    
    e = sub2ind([NELX_L, NELZ_B],  repmat([1:NELX_L]',iFZW,1), kron(NELZ_B - [0: iFZW - 1]', ones(NELX_L,1)));
    domain(1).etag(e)=tag_fz;
    
    e = sub2ind([NELX_M, NELZ_B],  repmat([1:NELX_M]',iFZW,1), kron(NELZ_B - [0: iFZW - 1]', ones(NELX_M,1)));
    domain(2).etag(e)=tag_fz;
    
    e = sub2ind([NELX_R, NELZ_B],  repmat([1:NELX_R]',iFZW,1), kron(NELZ_B - [0: iFZW - 1]', ones(NELX_R,1)));
    domain(3).etag(e)=tag_fz;
    
    e = sub2ind([NELX_R, NELZ_B],  repmat([1:NELX_R]',iFZW,1), kron([1: iFZW]', ones(NELX_R,1)));
    domain(4).etag(e)=tag_fz;
    
    e = sub2ind([NELX_M, NELZ_B],  repmat([1:NELX_M]',iFZW,1), kron([1: iFZW]', ones(NELX_M,1)));
    domain(5).etag(e)=tag_fz;
    
    e = sub2ind([NELX_L, NELZ_B],  repmat([1:NELX_L]',iFZW,1), kron([1: iFZW]', ones(NELX_L,1)));
    domain(6).etag(e)=tag_fz;
   
end
    
% tag the elements close to the fault
% to create a Kelvin-Voigt viscous layer
e = sub2ind([NELX_L, NELZ_B], NELX_L, NELZ_B);

if KVW>=h
    ikvw = ceil(KVW/h);
    
    [indx, indz] = meshgrid(NELX_L - [0 : ikvw - 1], NELZ_B - [0: ikvw - 1]);
    e = sub2ind([NELX_L, NELZ_B],  indx(:), indz(:));
    domain(1).etag(e)=tag_kv;
    
    [indx, indz] = meshgrid(1: NELX_M, NELZ_B - [0: ikvw - 1]);
    e = sub2ind([NELX_M, NELZ_B],  indx(:), indz(:));
    domain(2).etag(e)=tag_kv;
    
    [indx, indz] = meshgrid([1: ikvw], NELZ_B - [0: ikvw - 1]);
    e = sub2ind([NELX_R, NELZ_B],  indx(:), indz(:));
    domain(3).etag(e)=tag_kv;
    
    [indx, indz] = meshgrid([1: ikvw], [1: ikvw]);
    e = sub2ind([NELX_R, NELZ_T],  indx(:), indz(:));
    domain(4).etag(e)=tag_kv;
        
    [indx, indz] = meshgrid([1: NELX_M], [1: ikvw]);
    e = sub2ind([NELX_M, NELZ_T],  indx(:), indz(:));
    domain(5).etag(e)=tag_kv;
    
    [indx, indz] = meshgrid(NELX_L - [0 : ikvw - 1], [1: ikvw]);
    e = sub2ind([NELX_L, NELZ_T],  indx(:), indz(:));
    domain(6).etag(e)=tag_kv;

end

% merging domain into a single mesh.

% create a mesh tab
mtab = zeros(4, 6);
mtab(1, 1) =  -1;
mtab(1, 2) =  -1;
mtab(1, 3) =  -1;
mtab(1, 4) =   3;
mtab(1, 5) =  -6;
mtab(1, 6) =   1;

mtab(2, 1) =  2;
mtab(2, 2) =  3;
mtab(2, 3) = -2;
mtab(2, 4) = -2;
mtab(2, 5) =  4;
mtab(2, 6) =  5;

mtab(3, 1) =  6;
mtab(3, 2) = -5;
mtab(3, 3) =  4;
mtab(3, 4) = -3;
mtab(3, 5) = -3;
mtab(3, 6) = -3;

mtab(4, 1) =  -4;
mtab(4, 2) =  1;
mtab(4, 3) =  2;
mtab(4, 4) =  5;
mtab(4, 5) =  6;
mtab(4, 6) = -4;


% force merging two fault tip nodes

%               nod_mtab(:,4)   node merging table (optional):
%                               merge node #nod_mtab(i,2) of domain #nod_mtab(i,1)
%                               with  node #nod_mtab(i,4) of domain #nod_mtab(i,3)
node_mtab = zeros(2, 4);

% merge the left fault tip
d2 = domain(2);
d5 = domain(5);
d2_e_left  = sub2ind([NELX_M, NELZ_B],  1, NELZ_B);
tip_left_2 = d2.enod(4,d2_e_left);

d2_e_right  = sub2ind([NELX_M, NELZ_B], NELX_M, NELZ_B);
tip_right_2 = d2.enod(3,d2_e_right);

d5_e_left  = sub2ind([NELX_M, NELZ_T],  1, 1);
tip_left_5 = d5.enod(1,d5_e_left);

d5_e_right  = sub2ind([NELX_M, NELZ_T], NELX_M, 1);
tip_right_5 = d5.enod(2,d5_e_right);

node_mtab(1, :) = [2, tip_left_2, 5, tip_left_5];
node_mtab(2, :) = [2, tip_right_2, 5, tip_right_5];

mesh = mesh2d_merge(domain,mtab,node_mtab);

figure(2);
clf
mesh2d_plot(mesh, 2);
axis equal

% -- export the mesh file
mesh2d_write(mesh,'rsf_twosides');





