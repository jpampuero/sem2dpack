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
