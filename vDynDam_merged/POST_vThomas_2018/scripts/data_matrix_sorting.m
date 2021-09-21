% Script to rearrange the data in 2D matrix

% Marion Thomas last modified september 2019

%CALLS: 

%==========================================================================
%% DATA REARRANGEMENT

%Number of element:
Mx=max(grid.coord(:,1));
mx=min(grid.coord(:,1));
Mz=max(grid.coord(:,2));
mz=min(grid.coord(:,2));
Nz=abs(Mz-mz)/h;
Nx=abs(Mx-mx)/h;

% Data
target_box=[];
if ndims(v) == 3
    [wrk.indx,wrk.coord] = initialize(grid.ibool,grid.coord,target_box);
else
    wrk.indx = initialize(grid.ibool);
end
ll=size(wrk.indx(:,:),1);
cc=size(wrk.indx(:,:),2);
if ndims(v) == 3
    x=reshape(wrk.coord(wrk.indx(:,:),1),ll,cc)';
    z=reshape(wrk.coord(wrk.indx(:,:),2),ll,cc)';
else
    x=reshape(grid.coord(wrk.indx(:,:),1),ll,cc)';
    z=reshape(grid.coord(wrk.indx(:,:),2),ll,cc)';
end
field=reshape(v(wrk.indx(:,:)),ll,cc)';


%Reshape for pcolor and interpolation
Nline=(Mz-mz)/h;
Ncol=(Mx-mx)/h;
xrsh=[];zrsh=[];crsh=[];
for r=1:Nline/2
	xrsh=[xrsh,x(:,(r-1)*Ncol*2+1:2:r*Ncol*2),x(:,(r-1)*Ncol*2+2:2:r*Ncol*2)];
	zrsh=[zrsh,z(:,(r-1)*Ncol*2+1:2:r*Ncol*2),z(:,(r-1)*Ncol*2+2:2:r*Ncol*2)];
	crsh=[crsh,field(:,(r-1)*Ncol*2+1:2:r*Ncol*2),field(:,(r-1)*Ncol*2+2:2:r*Ncol*2)];
end 
xmat=reshape(mean(xrsh)',Ncol,Nline)';
zmat=reshape(mean(zrsh)',Ncol,Nline)';
cmat=reshape(mean(crsh)',Ncol,Nline)';

%% Resampling

gridsize = 0.05 * proZ; %1.0*R0;
[xsamp,zsamp] = meshgrid(mx:gridsize:Mx, mz:gridsize:Mz);
csamp = griddata(xmat,zmat,cmat,xsamp,zsamp, 'cubic');
Maxcsamp=max(max(csamp));
K = 0.15*ones(3);%K = 0.25*ones(3);
cconv = conv2(csamp./Maxcsamp,K,'same');
csamp=(cconv./(max(max(cconv))))*Maxcsamp;


