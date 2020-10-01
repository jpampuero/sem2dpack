% analyze test results for rsf with symmetry, same as Kaneko et al., 2011
% simulate half of the domain with absorbing bc for all external bcs.
close all; clear all;
addpath(genpath('../../../POST/'));
addpath(genpath('../../../PRE/'));

% read grid and plot grid
wdir = pwd;
grid = sem2d_read_specgrid(wdir);
h=sem2d_plot_grid(grid);
df = sem2d_read_fault('Flt05');
%% plot the slip distribution

figure(1);
plot(df.x, df.v(:, 1:3));
legend({'1','2','3'});

%% plot the slip velocity in the domain
it=30;
vs = sem2d_snapshot_read('vy', it, wdir);
figure(2)
sem2d_snapshot_plot(vs, grid);
% set(gca,'DefaultPatchEdgeColor','k');
hold on;
h=sem2d_plot_grid(grid);
hold off;
% ylim([0, 5])
title('v');
set(gca,'fontsize',14);