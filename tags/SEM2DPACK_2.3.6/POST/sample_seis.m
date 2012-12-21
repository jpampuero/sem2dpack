% SAMPLE_SEIS example of visualization of seismogram output from SEM2DPACK
% Assumes a 2D P-SV simulation (ndof=2)

% Read seismogram outputs from SEM2DPACK
data = sem2d_read_seis();

% Plot all x-component traces together, offset by distance to first station
figure(1)
d = sqrt((data.x-data.x(1)).^2 +(data.z-data.z(1)).^2);
plot_seis(d,data.dt,data.ux)

% Plot seismograms (two components) at first station
t = [1:data.nt]*data.dt;
figure(2)
subplot(211)
plot(t,data.ux(:,1))
ylabel('Ux')
subplot(212)
plot(t,data.uz(:,1))
ylabel('Uz')
xlabel('Time (s)')
