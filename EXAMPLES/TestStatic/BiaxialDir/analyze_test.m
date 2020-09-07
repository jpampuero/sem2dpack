% plot biaxial dirichlet loading

% read grid and plot grid
wdir = pwd;
grid = sem2d_read_specgrid(wdir);
% h=sem2d_plot_grid(grid);

%% analytical solution
% e11 = 0, e22 = 1e-3, e33 = 0, K = 3, Mu = 1

vp = 1.7321; rho=1; vs = 1;

% elastic properties
Mu = vs^2*rho;
K = vp^2*rho - 4/3*Mu;

E = 9*K*Mu/(3*K + Mu);
Lam =  K - 2*Mu/3;
% Nu = (3*K - 2*Mu)/2/(3*K + Mu);

% Lam=1;
e_11 = -2e-3;
e_22 = 1e-3;

s11a = (Lam + 2 * Mu) * e_11 + Lam * e_22;
s22a = Lam * e_11 + (Lam +  2 * Mu) * e_22;

%%
it = 1;

% read fields
% dx = sem2d_snapshot_read('dx', it, wdir);
% dz = sem2d_snapshot_read('dz', it , wdir);

% vx = sem2d_snapshot_read('vx', it, wdir);
% vz = sem2d_snapshot_read('vz', it , wdir);

% read stresses and strains
% e11 = sem2d_snapshot_read('e11', it , wdir);
s11 = sem2d_snapshot_read('s11', it , wdir);

% e22 = sem2d_snapshot_read('e22', it , wdir);
s22 = sem2d_snapshot_read('s22', it , wdir);

% sem2d_snapshot_plot(e11, grid, [0 1e-3]);
assert(abs(mean(s11(:)) - s11a)/s11a<1e-6);
assert(abs(mean(s22(:)) - s22a)/s22a<1e-6);

fprintf('Test = 1\n');
fprintf('Biaxial Loading: TEST PASSED!\n')
fprintf('s11 sem: %8.6f, s11 analytical: %8.6f\n', mean(s11(:)) , s11a);
fprintf('s22 sem: %8.6f, s22 analytical: %8.6f\n', mean(s22(:)) , s22a);
