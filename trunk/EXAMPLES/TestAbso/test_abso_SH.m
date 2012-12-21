

% Read parameters from header file
[dt,nsamp,nsta] = textread('SeisHeader_sem2d.hdr','%n%n%n',1,'headerlines',1);
[xsta,zsta] = textread('SeisHeader_sem2d.hdr','%f%f','headerlines',3);

% Read seismograms
fid=fopen('Uy1_sem2d.dat'); uy1 = fread(fid,[nsamp,nsta],'single') ; fclose(fid);
fid=fopen('Uy2_sem2d.dat'); uy2 = fread(fid,[nsamp,nsta],'single') ; fclose(fid);
%uy1=cumsum(uy1); uy2=cumsum(uy2);

t = [1:nsamp]*dt;

% offset by x coordinate
doff = abs(xsta(end)-xsta(1))/(nsta-1);
offset = repmat( xsta', nsamp,1); 

a_factor = 3; % factor for seismogram amplitude (->overlap)
m_factor = 10; % factor for misfit

t0 = 0.05; % source time shift
delta = 1e3; % distance from source to receiver line
vp = 5000; %4000; % P wave velocity
vs = 2887; %2310; % S wave velocity

%----
if PLOT_FIGURE_1

% Plot all Uy traces together
ascale = a_factor*doff/max(abs(uy2(:)));
figure(1)
clf
subplot(121)
plot(t,uy1*ascale +offset,'--', t,uy2*ascale +offset)
xlabel('Time (s)')
ylabel('Uy -- Distance (m)')
title('Seismograms')
subplot(122)
plot(t,(uy2-uy1)*ascale*m_factor +offset)
xlabel('Time (s)')
title('Error *10')

end

%----
if PLOT_FIGURE_2

figure(2)
clf
plot(t0+sqrt(delta^2+xsta.^2)/vs,xsta,'.-', ... % S from source
     t0+sqrt((5*delta)^2+xsta.^2)/vs,xsta,'-+', ... % S from top & bottom
     t0+sqrt(delta^2+(3*delta+xsta).^2)/vs,xsta,'--', ... % S from left
     t0+sqrt(delta^2+(-9*delta+xsta).^2)/vs,xsta,'--' ) % S from right
hold on
ascale = a_factor*doff/max(abs(uy2(:)));
plot(t,uy2*ascale*a_factor +offset)
hold off
axis([0 2 -500 5000])
legend('S from source','S from top & bottom','S from left','S from right',2)

end

%---

if PLOT_FIGURE_3

figure(3)
clf

%-- read previous results
%load stacey_old % theta errp errs
%load paraxial_explo % theta errp errs
%plot(theta,errs)

theta = atan(xsta/delta)*180/pi;
tw = 1.5*t0;

% analyze direct S
tp = t0+sqrt(delta^2+xsta.^2)/vs;
win1 = floor( (tp-tw)/dt );
win2 = ceil( (tp+tw)/dt );
errs=zeros(nsta,1);
for ista=1:nsta,
  if win2(ista)>nsamp, break, end
  win = [win1(ista):win2(ista)];
  errs(ista) = sqrt( ...
    sum((uy2(win,ista)-uy1(win,ista)).^2) ./ sum(uy2(win,ista).^2)  );
end

hold on
plot(theta, errs,'o', theta, (1-cos(pi*theta/180))./(1+cos(pi*theta/180)) , 'k--')
xlabel('\theta (degrees)')
ylabel('Relative misfit')
hold off

end
