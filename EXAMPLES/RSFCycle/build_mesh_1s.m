% build mesh for rate state cycle simulation
%
% One-sided faults
%
% a finite single fault in the middle of the domain
%
%   | --------------------------------------------------|
%   |                                                               |
%   |                                                               |
%   |                     ----------------                      |
%
clear; close all;
addpath(genpath('../../PRE/'));

%
Lx    =   91;
Lz    =   60;

%-- you can modify the following parameters:
FAULT_L    = 45; % fault length 	
FAULT_Z    = 0;   % fault z location 
FZW           = 0;   % fault weak zone thickness each side
KVW           = 0;  % width of KV material on each side
tag_fz         = 2;   % tag for fault zone
tag_kv        = 3;   % tag for KV damper, tag_kv=0, no KV damper
h                = 1;    % element size

% ----------------------------- Do not modify --------------------------
% create vertices
p1   = [-Lx/2,              0];
p2   = [-FAULT_L/2,   0];
p3   = [FAULT_L/2,    0];
p4   = [ Lx/2,               0];

p5   = [-Lx/2,              Lz];
p6  = [-FAULT_L/2,   Lz];
p7  = [FAULT_L/2,    Lz];
p8  = [Lx/2,               Lz];

% number of elements for left, middle, right columns
NELX_L  =  ceil((Lx - FAULT_L)/2/h);
NELX_M = ceil(FAULT_L/h);
NELX_R =  NELX_L;

NELZ_T = ceil(Lz/h);

% create  domains

% LEFT
ps = [p1; p2; p6; p5];
domain(1) = mesh2d_quad(ps(:,1), ps(:,2), NELX_L, NELZ_T);

% MID
ps = [p2; p3; p7; p6];
domain(2) = mesh2d_quad(ps(:,1), ps(:,2), NELX_M, NELZ_T);

% RIGHT
ps = [p3; p4; p8; p7];
domain(3) = mesh2d_quad(ps(:,1), ps(:,2), NELX_R, NELZ_T);

% tag the fault zone
if FZW >= h
    iFZW = ceil(FZW/h);
    
    e = sub2ind([NELX_R, NELZ_T],  repmat([1:NELX_R]',iFZW,1), kron([1: iFZW]', ones(NELX_R,1)));
    domain(1).etag(e)=tag_fz;
    
    e = sub2ind([NELX_M, NELZ_T],  repmat([1:NELX_M]',iFZW,1), kron([1: iFZW]', ones(NELX_M,1)));
    domain(2).etag(e)=tag_fz;
    
    e = sub2ind([NELX_L, NELZ_T],  repmat([1:NELX_L]',iFZW,1), kron([1: iFZW]', ones(NELX_L,1)));
    domain(3).etag(e)=tag_fz;
   
end
       23  iterations.
 quasi-static solve, pass            2
 PCG solver converges in          449  iterations.
Timestep #      30  t =   29.457E-03  vmax =  945.540E-06  dmax =   20.680E-06
 dt =    5.9163567141755297E-003
 quasi-static solve, pass            1

% tag the elements close to the fault
% to create a Kelvin-Voigt viscous layer

if KVW>=h
    ikvw = ceil(KVW/h);
    
    [indx, indz] = meshgrid([1: ikvw], [1: ikvw]);
    e = sub2ind([NELX_R, NELZ_T],  indx(:), indz(:));
    domain(1).etag(e)=tag_kv;
        
    [indx, indz] = meshgrid([1: NELX_M], [1: ikvw]);
    e = sub2ind([NELX_M, NELZ_T],  indx(:), indz(:));
    domain(2).etag(e)=tag_kv;
    
    [indx, indz] = meshgrid(NELX_L   23  iterations.
 quasi-static solve, pass            2
 PCG solver converges in          449  iterations.
Timestep #      30  t =   29.457E-03  vmax =  945.540E-06  dmax =   20.680E-06
 dt =    5.9163567141755297E-003
 quasi-static solve, pass            1
 - [0 : ikvw - 1], [1: ikvw]);
    e = sub2ind([NELX_L, NELZ_T],  indx(:), indz(:));
    domain(3).etag(e)=tag_kv;

end

% plot the unmerged domain
figure(1);
mesh2d_plot(domain, 2);
axis equal
%% merging domain into a single mesh.

% create a mesh tab
mtab = zeros(4, 3);
mtab(1, 1) =  -1;
mtab(1, 2) =  -5;
mtab(1, 3) =  -1;

mtab(2, 1) =  2;
mtab(2, 2) =  3;
mtab(2, 3) = -2;

mtab(3, 1) =  -3;
mtab(3, 2) =  -3;
mtab(3, 3) =  -3;

mtab(4, 1) = -4;
mtab(4, 2) =  1;
mtab(4, 3) =  2;

mesh = mesh2d_merge(domain,mtab);
figure(2);
clf
mesh2d_plot(mesh, 2);
axis equal

% -- export the mesh file
% mesh2d_write(mesh,'rsf_oneside');





