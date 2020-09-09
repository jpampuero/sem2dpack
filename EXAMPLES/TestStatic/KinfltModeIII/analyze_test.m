% plot kinematic fault test mode 3

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
it = 4;
fault_z = 0;

% obtain fault shear stress change from stresses
dx = sem2d_snapshot_read('dy', it, wdir);
s11 = sem2d_snapshot_read('s23', it , wdir);

%%
figure(1)
subplot(1, 2, 1);
sem2d_snapshot_plot(dx, grid, [-2, 2]*1e-4);
title('uy');
set(gca,'fontsize',14)
subplot(1,2, 2);
wrk = sem2d_snapshot_plot(s11, grid, [-2, 2]*1e-4);
title('syz');
set(gca,'fontsize',14)
% subplot(1,2, 3);
% wrk = sem2d_snapshot_plot(s22, grid, [-2, 2]*1e-4);
% title('szz');
cmap();

% plot fault stress change.
figure(2)
subplot(2,1,1);
plot(d_flt.x, d_flt.d(:, 5),'linew',1.5);
set(gca,'fontsize',16)
xlabel('x');
ylabel('slip');

s11 = s11(:);
flt_ind   = find(abs(wrk.coord(:,2))<1e-6);
xz_flt    = wrk.coord(flt_ind, :);
s11_flt  = s11(flt_ind);

% obtain unique x index
[x,IA,~] = unique(xz_flt(:, 1));
s11_flt = s11_flt(IA);
subplot(2,1,2);
plot(x, s11_flt,'linew',1.5);
set(gca,'fontsize',16)
xlabel('x');
ylabel('stress change');
