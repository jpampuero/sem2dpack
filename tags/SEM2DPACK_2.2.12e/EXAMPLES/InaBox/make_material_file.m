% generates a heterogeneous rho, vp ,vs model in a regular grid
% to be used in the DIST_HETE1 blocks

ncol=3;
nx=31;
nz=31;
x0=-1e3;
z0=-1e3;
dx=2e3/(nx-1);
dz=2e3/(nz-1);

% heterogeneous material model: background + random
rho = 2400*ones(1,nx*nz);
vp  = 3800*(1+0.3*(rand(1,nx*nz)-1/2)*2);
vs  = 2000/3800*vp;

fid=fopen('material.tab','w');
fprintf(fid,'%u %u %u %f %f %f %f\n',ncol,nx,nz,x0,z0,dx,dz);
fprintf(fid,'%f %f %f\n',[rho;vp;vs]);
fclose(fid);
