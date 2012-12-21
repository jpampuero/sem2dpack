PLOT_FIGURE_1 = 1;
PLOT_FIGURE_2 = 1;
PLOT_FIGURE_3 = 1;

% Read parameters from header file
[dt,nsamp,nsta] = textread('SeisHeader_sem2d.hdr','%n%n%n',1,'headerlines',1);
[xsta,zsta] = textread('SeisHeader_sem2d.hdr','%f%f','headerlines',3);

% Read seismograms
fid=fopen('Ux1_sem2d.dat'); ux1 = fread(fid,[nsamp,nsta],'single') ; fclose(fid);
fid=fopen('Uz1_sem2d.dat'); uz1 = fread(fid,[nsamp,nsta],'single') ; fclose(fid);
fid=fopen('Ux2_sem2d.dat'); ux2 = fread(fid,[nsamp,nsta],'single') ; fclose(fid);
fid=fopen('Uz2_sem2d.dat'); uz2 = fread(fid,[nsamp,nsta],'single') ; fclose(fid);

% if displacement instead of velocity:
ux1=cumsum(ux1); uz1=cumsum(uz1); ux2=cumsum(ux2); uz2=cumsum(uz2);

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

% Plot all Ux traces together
ascale = a_factor*doff/max(abs(ux2(:)));
figure(1)
clf
subplot(221)
plot(t,ux1*ascale +offset,'--', t,ux2*ascale +offset)
ylabel('Ux -- Distance (m)')
title('Seismograms')
subplot(222)
plot(t,(ux2-ux1)*ascale*m_factor +offset)
title('Error *10')

% Plot all Uz traces together
ascale = a_factor*doff/max(abs(uz2(:)));
subplot(223)
plot(t,uz1*ascale +offset,'--', t,uz2*ascale +offset)
xlabel('Time (s)')
ylabel('Uz -- Distance (m)')
subplot(224)
plot(t,(uz2-uz1)*ascale*m_factor +offset)
xlabel('Time (s)')

end

%----
if PLOT_FIGURE_2

figure(2)
clf
plot(t0+sqrt(delta^2+xsta.^2)/vp,xsta,'.-', ... % P from source
     t0+sqrt((5*delta)^2+xsta.^2)/vp,xsta,'-+', ... % P from top & bottom
     t0+sqrt(delta^2+(3*delta+xsta).^2)/vp,xsta,'--', ... % P from left
     t0+sqrt(delta^2+(-9*delta+xsta).^2)/vp,xsta,'--', ... % P from right
     t0+sqrt(delta^2+xsta.^2)/vs,xsta,'.-') % S from source
%     t0+sqrt((5*delta)^2+xsta.^2)/vs,xsta,'.-') % S from top & bottom
hold on
ip = [16:24]*pi/180;
is = asin(vs/vp *sin(ip));
x = (delta -1.5*delta*tan(ip)) ./tan(is)  ; % |x| from boundary
plot( t0+ 1.5*delta/vp./cos(ip) + x/vs./cos(is), x - 1.5*delta, '-x')
ip = [7.5:12.5]*pi/180;
is = asin(vs/vp *sin(ip));
x = (delta -4.5*delta*tan(ip)) ./tan(is)  ; % |x| from boundary
plot( t0+ 4.5*delta/vp./cos(ip) + x/vs./cos(is), -x + 4.5*delta, 'g-x')
%plot(t,sqrt(ux2.^2+uz2.^2)*ascale +offset)
ascale = a_factor*doff/max(abs(ux2(:)+uz2(:)));
plot(t,(ux2+uz2)*ascale*a_factor +offset)
%plot(t,ux2*ascale +offset,t,uz2*ascale +offset)
hold off
axis([0 1.5 -500 5000])
legend('P from source','P from top & bottom','P from left','P from right',...
       'S from source', ...
       'P-S from left','P-S from right',2)

end

%---

if PLOT_FIGURE_3

figure(3)
clf

%-- read previous results
%load stacey_old % theta errp errs
load paraxial_explo % theta errp errs
plot(theta,errp,theta(1:length(errs)),errs)

theta = atan(xsta/delta)*180/pi;
tw = 1.5*t0;

% analyze direct P
tp = t0+sqrt(delta^2+xsta.^2)/vp;
win1 = floor( (tp-tw)/dt );
win2 = ceil( (tp+tw)/dt );
errp=zeros(nsta,1);
for ista=1:nsta,
  if win2(ista)>nsamp, break, end
  win = [win1(ista):win2(ista)];
  errp(ista) = sqrt( ...
    sum( (ux2(win,ista)-ux1(win,ista)).^2 + (uz2(win,ista)-uz1(win,ista)).^2 ) ...
    ./ sum( ux2(win,ista).^2 + uz2(win,ista).^2 )  );
end

% analyze direct S
tp = t0+sqrt(delta^2+xsta.^2)/vs;
win1 = floor( (tp-tw)/dt );
win2 = ceil( (tp+tw)/dt );
% keep only direct S arriving well before the P reflected from top & bottom
tp2 = t0+sqrt((5*delta)^2+xsta.^2)/vp; tp2=floor((tp2-tw)/dt); 
nstaS = nnz( win2 <= tp2 );
errs=zeros(nstaS,1);
for ista=1:nstaS,
  if win2(ista)>nsamp, break, end
  win = [win1(ista):win2(ista)];
  errs(ista) = sqrt( ...
    sum( (ux2(win,ista)-ux1(win,ista)).^2 + (uz2(win,ista)-uz1(win,ista)).^2 ) ...
    ./ sum( ux2(win,ista).^2 + uz2(win,ista).^2 )  );
end

hold on
plot(theta,errp,'.-', theta(1:nstaS), errs,'.-')
xlabel('\theta (degrees)')
ylabel('Relative misfit')
legend('P','S',2)
hold off

end
