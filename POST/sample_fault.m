% Read parameters from header file
[npts,ndat,nsamp,dt] = textread('Flt01_sem2d.hdr','%n%n%n%n',1,'headerlines',1);
[x,z] = textread('Flt01_sem2d.hdr','%f%f','headerlines',4);

% Read fault data in a big matrix
fid=fopen('Flt01_sem2d.dat'); raw = fread(fid,[npts+2,inf],'single') ; fclose(fid);
raw = reshape(raw(2:npts+1,:),[npts ndat nsamp]);

% reformat each field [npts,nsamp]
Slip = squeeze(raw(:,1,:)); 
SlipRate = squeeze(raw(:,2,:)); 
ShearStress = squeeze(raw(:,3,:)); 
NormalStress = squeeze(raw(:,4,:)); 
Friction = squeeze(raw(:,5,:)); 

% Plot slip rate
mesh([0:nsamp-1]*dt,x/1e3,SlipRate)
xlabel('Time (s)')
ylabel('Along strike distance (km)')
zlabel('Slip rate (m/s)')
