P_OR_S = 1;
T0 = 0.05;
TW = T0;
VP = 5e3;
VS = 2887;
L = 500;

% Read parameters from header file
[dt,nsamp,nsta] = textread('SeisHeader_sem2d.hdr','%n%n%n',1,'headerlines',1);
[xsta,zsta] = textread('SeisHeader_sem2d.hdr','%f%f','headerlines',3);

% Read seismograms
fid=fopen('Ux_sem2d.dat'); ux = fread(fid,[nsamp,nsta],'single') ; fclose(fid);
fid=fopen('Uz_sem2d.dat'); uz = fread(fid,[nsamp,nsta],'single') ; fclose(fid);

t = [1:nsamp]*dt;

% offset by x coordinate
doff = abs(xsta(end)-xsta(1))/(nsta-1);
offset = repmat( xsta', nsamp,1); 

if P_OR_S==1, 
  V=VP; u=ux; 
else
  V=VS; u=uz; 
end

% Plot all traces together
figure(1)
clf
ascale = doff/max(abs(u(:)));
plot(t,u*ascale +offset, T0+xsta/V,xsta,'.--',T0 - (xsta-2*L)/V,xsta,'.--' )
xlabel('Time (s)')
ylabel('Distance (m)')

% Plot error / frequency
p1 = floor( (T0 + xsta/V - TW)/dt )+1;
p2 = floor( (T0 - (xsta-2*L)/V -TW)/dt )+1;
tw = 2*floor(TW/dt);
ista = 1;
nfft = 2^nextpow2(tw);
ferr = abs( fft(u(p2(ista):p2(ista)+tw,ista),nfft) ) ...
     ./ abs( fft(u(p1(ista):p1(ista)+tw,ista),nfft) );
ferr = ferr(1:nfft/2+1);
f = [0:nfft/2] /(nfft*dt);
figure(2)
clf
loglog(f,ferr)
err = norm(u(p2(ista):p2(ista)+tw,ista)) / norm(u(p1(ista):p1(ista)+tw,ista));
title( sprintf('Relative RMS error = %g',err) )
xlabel('Frequency (Hz)')
ylabel('Relative error')
grid on; grid minor
