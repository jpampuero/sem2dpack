% Read parameters from header file
[dt,nsamp,nsta] = textread('SeisHeader_sem2d.hdr','%n%n%n',1,'headerlines',1);
[xsta,zsta] = textread('SeisHeader_sem2d.hdr','%f%f','headerlines',3);

% Read seismograms
fid=fopen('Ux_sem2d.dat'); ux = fread(fid,[nsamp,nsta],'single') ; fclose(fid);
fid=fopen('Uz_sem2d.dat'); uz = fread(fid,[nsamp,nsta],'single') ; fclose(fid);

t = [1:nsamp]*dt;

% Plot all Ux traces together, offset by distance to first station
figure(1)
d = sqrt((xsta-xsta(1)).^2 +(zsta-zsta(1)).^2);
doff = max(d)/(nsta-1);
ascale = doff/max(abs(ux(:)));
offset = [0:nsta-1]*doff; 
plot(t,ux*ascale +repmat(offset,nsamp,1))
xlabel('Time (s)')
ylabel('Distance (m)')

% Plot fisrt station
figure(2)
subplot(212)
plot(t,ux(:,1))
ylabel('Ux')
subplot(212)
plot(t,uz(:,1))
ylabel('Uz')
xlabel('Time (s)')
