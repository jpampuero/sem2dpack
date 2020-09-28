% plot kinematic bc test
addpath(genpath('../../../POST/'));
addpath(genpath('../../../PRE/'));

% read grid and plot grid
wdir = pwd;
grid = sem2d_read_specgrid(wdir);
h=sem2d_plot_grid(grid);
d_flt = sem2d_read_kinflt('Flt05');

%% analytical solution
vp = 1.7321; rho=1; vs = 1;

% elastic properties
Mu = vs^2*rho;
K = vp^2*rho - 4/3*Mu;

E = 9*K*Mu/(3*K + Mu);
Lam =  K - 2*Mu/3;
Nu = (3*K - 2*Mu)/2/(3*K + Mu);

%% plot fault slip and stress change
close all;
it = 1;
fault_z = 0;

% obtain fault shear stress change from stresses
vy = sem2d_snapshot_read('vy', it, wdir);

figure(1)
sem2d_snapshot_plot(vy, grid);
title('vy');
set(gca,'fontsize',14)