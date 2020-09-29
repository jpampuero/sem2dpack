% analyze test for Agelaw with constant stress.
% results should produce linear increase in state variable theta
% and hyperbolic decay of slip velocity over time
% 
% Due to the boundary condition, rigid body motion exists in the solution.
%
close all;
addpath(genpath('../../../POST/'));
addpath(genpath('../../../PRE/'));

% read grid and plot grid
wdir = pwd;
grid = sem2d_read_specgrid(wdir);
h=sem2d_plot_grid(grid);
d_flt = sem2d_read_fault('Flt05');
mode  = 2;

%% compute analytical results

a   = 0.01;
f0  = 0.6;
V0 = 1e-6;
f    = 0.55; % f<f0, v increases with time.
L   = 1e-4;
C   = exp((f - f0)/a);

% initial velocity
V_i = 1e-3;
theta_0 = C*L/V_i;

t    = 0:0.0001:0.0020;
theta = (1 - C) * t + theta_0;
V    = C*L./theta;
%% compare analytical and numerical results
V_n  = mean(d_flt.v, 1);
theta_n = mean(d_flt.theta, 1);

figure(1);
subplot(1,2,1);
plot(V,'-','linew',1.5);
hold on;
plot(V_n, '*');
hold off;
xlabel('step');
ylabel('V');
set(gca,'fontsize',14);
legend({'analytic','sem'});
title('slip velocity');

subplot(1,2,2);
plot(theta,'-','linew',1.5);
hold on;
plot(theta_n, '*');
hold off;
xlabel('step');
ylabel('theta');
set(gca,'fontsize',14);
legend({'analytic','sem'});
title('state');

%% plot the velocity field 
it = 19;

% obtain fault shear stress change from stresses
switch mode
    case 3
    vs = sem2d_snapshot_read('vy', it, wdir);
    vs = vs-mean(vs(:)); % remove the rigid body motion
    case 2
    vs = sem2d_snapshot_read('vx', it, wdir);
    vs = vs-mean(vs(:)); % remove the rigid body motion
end

figure(3)
sem2d_snapshot_plot(vs, grid);
set(gca,'DefaultPatchEdgeColor','k');
hold on;
h=sem2d_plot_grid(grid);
hold off;
title('v - mean');
set(gca,'fontsize',14);shg

fprintf('Step %d, slip v from vs: %f\n', it, max(vs(:))-min(vs(:)) );
fprintf('Step %d, slip v saved on fault: %f\n', it-1, mean(d_flt.v(:, it-1)));
fprintf('Step %d, slip v saved on fault: %f\n', it, mean(d_flt.v(:, it)));
fprintf('Step %d, slip v saved on fault: %f\n', it+1, mean(d_flt.v(:, it+1)));

