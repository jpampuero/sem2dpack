% analyze test results for rsf with symmetry, same as Kaneko et al., 2011
% simulate half of the domain with absorbing bc for all external bcs.
close all; clear all;
addpath(genpath('../../../POST/'));
addpath(genpath('../../../PRE/'));

% read grid and plot grid
wdir = pwd;
grid = sem2d_read_specgrid(wdir);
df = sem2d_read_fault('Flt05',wdir);

wdir = '../kaneko11_homo_sym';
grid1 = sem2d_read_specgrid(wdir);
df1 = sem2d_read_fault('Flt05',wdir);
%% plot the slip distribution

figure(1);
plot(df1.x, df1.v(:, 1:15),'k-');
hold on;
plot(df.x, df.v(:, 1:15),'r-');
hold off;

%% plot the slip velocity in the domain
it=30;
vs = sem2d_snapshot_read('vy', it, wdir);
figure(2)
sem2d_snapshot_plot(vs, grid);
%% set(gca,'DefaultPatchEdgeColor','k');
%hold on;
%h=sem2d_plot_grid(grid);
%hold off;
%% ylim([0, 5])
%title('v');
%set(gca,'fontsize',14);
