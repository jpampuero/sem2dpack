% compare SEM2DPACK seismogram to analytic seismogram
% for SH point force source

% read SEM2DPACK seismograms
[dt,nsamp,nsta] = textread('SeisHeader_sem2d.hdr','%n%n%n',1,'headerlines',1);
t = [0:nsamp-1]*dt;
fid=fopen('Uy_sem2d.dat'); u = fread(fid,[nsamp,nsta],'single') ; fclose(fid);
u = u(:,5);

% read analytical seismograms
% generated with analytic_sh2d.m (SEMLAB)
load uyref

err = max(abs( u-uref )) ./ max(abs( uref ));
test = max(err)<0.02;

disp(['L1 misfit uy = ' num2str(err)])
disp(['Test = ' num2str(test) ])

%plot(t,uy, tref,uyref, t,(uy-uyref)*10)
