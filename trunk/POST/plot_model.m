% PLOT_MODEL plots velocity and density model

grid = sem2d_read_specgrid();

subplot(311)
cp = sem2d_snapshot_read('Cp');
sem2d_snapshot_plot(cp,grid,[0 inf]);
axis equal
title('Cp (m/s)')

subplot(312)
cs = sem2d_snapshot_read('Cs');
sem2d_snapshot_plot(cs,grid,[0 inf]);
axis equal
title('Cs (m/s)')

subplot(313)
rho = sem2d_snapshot_read('Rho');
sem2d_snapshot_plot(rho,grid,[0 inf]);
axis equal
title('density (kg/m^3)')

