% SAMPLE_FAULT example of visualization of fault output from SEM2DPACK

% Read fault data
data = sem2d_read_fault('Flt01');

% Plot slip rate
mesh([0:data.nt-1]*data.dt,data.x/1e3,data.v)
xlabel('Time (s)')
ylabel('Along strike distance (km)')
zlabel('Slip rate (m/s)')
