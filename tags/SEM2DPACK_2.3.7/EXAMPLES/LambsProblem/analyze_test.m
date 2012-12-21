% compare SEM2DPACK seismogram to analytic seismogram
% for Lamb's problem

% read SEM2DPACK seismograms
[dt,nsamp,nsta] = textread('SeisHeader_sem2d.hdr','%n%n%n',1,'headerlines',1);
t = [0:nsamp-1]*dt;
fid=fopen('Ux_sem2d.dat'); ux = fread(fid,[nsamp,nsta],'single') ; fclose(fid);
fid=fopen('Uz_sem2d.dat'); uz = fread(fid,[nsamp,nsta],'single') ; fclose(fid);

% read analytical seismograms
% generated with EX2DDIR (from the SPICE software repository)
uxa = load('Ux_file_ascii');
uxa = reshape(uxa,length(uxa)/nsta,nsta);
uza = load('Uz_file_ascii');
uza = reshape(uza,length(uza)/nsta,nsta);
% ta = [1:size(uza,1)]*dt;
% add time=0
uxa = [zeros(1,nsta) ; uxa];
uza = [zeros(1,nsta) ; uza];

err = max(abs( [ux-uxa , uz-uza] )) ./ max(abs( [uxa , uza] ));
test = max(err)<0.005;

disp(['L1 misfit ux1 ux2 uz1 uz2 = ' num2str(err)])
disp(['Test = ' num2str(test) ])
